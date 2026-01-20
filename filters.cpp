#include "plugin.h"
#include <cmath>

// Cross-platform SIMD support
#ifdef __APPLE__
    #if defined(__x86_64__) || defined(__i386__)
        #include <immintrin.h>  // x86 SIMD
    #elif defined(__aarch64__) || defined(__arm64__)
        // ARM64: disable x86 SIMD for now
        #undef __SSE__
    #endif
#else
    #include <immintrin.h>
#endif

// Constant power panning function
void calc_constant_power_pan(float spread, float *pan_L, float *pan_R) {
    // Convert spread (-1 to 1) to angle (0 to PI/2)
    float angle = (spread + 1.0f) * (M_PI / 4.0f);  // 0 to PI/2
    
    *pan_L = cosf(angle);
    *pan_R = sinf(angle);
    
    // Apply minimum level for spread effect, but maintain constant power
    *pan_L = fmaxf(0.1f, *pan_L);
    *pan_R = fmaxf(0.1f, *pan_R);
    
    // Normalize to maintain constant total power
    float total_power = (*pan_L * *pan_L) + (*pan_R * *pan_R);
    float norm_factor = 1.0f / sqrtf(total_power);
    
    *pan_L *= norm_factor;
    *pan_R *= norm_factor;
}

// Optimized single filter processing with unrolled loop for common cases
static inline float process_single_biquad(float input, 
                                          float b0, float b1, float b2, float a1, float a2,
                                          float &x1, float &x2, float &y1, float &y2) {
    float output = b0 * input + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
    
    x2 = x1;
    x1 = input;
    y2 = y1;
    y1 = output;
    
    return output;
}

// Process a single channel through all active filters with optimization
void process_channel_filters(minimal_plugin_t *p, float *buffer, uint32_t frames, bool is_left_channel) {
    // Early exit if no active filters
    if (p->num_active_filters == 0) {
        return;
    }
    
    // Calculate constants once
    float hidden_boost = powf(10.0f, 12.0f / 20.0f);  // +12dB hidden boost
    float wet_boost_linear = hidden_boost * powf(10.0f, p->wet_boost_smooth.current / 20.0f);
    const float oversample_compensation = 1.0f;  // +0dB compensation
    
    float pan_factor_L, pan_factor_R;
    calc_constant_power_pan(p->spread_smooth.current, &pan_factor_L, &pan_factor_R);
    float pan_factor = is_left_channel ? pan_factor_L : pan_factor_R;
    
    // Get filter state pointers for this channel
    float (*bp_x1)[NUM_SEMITONES] = is_left_channel ? p->bp_x1_L : p->bp_x1_R;
    float (*bp_x2)[NUM_SEMITONES] = is_left_channel ? p->bp_x2_L : p->bp_x2_R;
    float (*bp_y1)[NUM_SEMITONES] = is_left_channel ? p->bp_y1_L : p->bp_y1_R;
    float (*bp_y2)[NUM_SEMITONES] = is_left_channel ? p->bp_y2_L : p->bp_y2_R;
    
    float &lpf_state = is_left_channel ? p->lpf_state_L : p->lpf_state_R;
    float &lpf2_state = is_left_channel ? p->lpf2_state_L : p->lpf2_state_R;
    float &lpf3_state = is_left_channel ? p->lpf3_state_L : p->lpf3_state_R;
    
    // Saturation parameters
    float internal_gain = 1.2f;
    float sat_amount = p->saturation_smooth.current * internal_gain;
    bool bypass_saturation = sat_amount < 0.001f;
    float inv_internal_gain = 1.0f / internal_gain;
    
    // Dry/wet mix parameters
    float dry_gain = (1.0f - p->dry_wet_smooth.current) * oversample_compensation;
    float wet_gain = p->dry_wet_smooth.current * wet_boost_linear * oversample_compensation * pan_factor;
    
    // Unroll filter processing for common cases to reduce loop overhead
    switch (p->num_active_filters) {
        case 1: {
            // Single filter - most common case, fully optimized
            int oct = p->active_filters[0].octave_idx;
            int semi = p->active_filters[0].semitone_idx;
            
            float b0 = p->bp_b0[oct][semi];
            float b1 = p->bp_b1[oct][semi];
            float b2 = p->bp_b2[oct][semi];
            float a1 = p->bp_a1[oct][semi];
            float a2 = p->bp_a2[oct][semi];
            
            for (uint32_t frame = 0; frame < frames; ++frame) {
                float input_sample = buffer[frame];
                
                float bp_output = process_single_biquad(input_sample, b0, b1, b2, a1, a2,
                                                       bp_x1[oct][semi], bp_x2[oct][semi], 
                                                       bp_y1[oct][semi], bp_y2[oct][semi]);
                
                // Apply saturation
                float sat_output = bypass_saturation ? bp_output : 
                                  tanhf(bp_output * sat_amount) * inv_internal_gain;
                
                // Apply LPF and mix
                float lpf_output = process_adaptive_lpf(sat_output, lpf_state, lpf2_state, lpf3_state, 
                                                       p->lpf_coeff, p->lpf_cutoff_smooth.current);
                
                buffer[frame] = input_sample * dry_gain + lpf_output * wet_gain;
                
                // Safety clipper
                if (fabsf(buffer[frame]) > 0.95f) {
                    buffer[frame] = tanhf(buffer[frame] * 0.7f) / 0.7f;
                }
            }
            break;
        }
        
        case 2: {
            // Two filters - second most common case
            for (uint32_t frame = 0; frame < frames; ++frame) {
                float input_sample = buffer[frame];
                float combined_output = 0.0f;
                
                // Unroll the two filter operations
                for (int i = 0; i < 2; i++) {
                    int oct = p->active_filters[i].octave_idx;
                    int semi = p->active_filters[i].semitone_idx;
                    
                    combined_output += process_single_biquad(input_sample,
                                                           p->bp_b0[oct][semi], p->bp_b1[oct][semi], p->bp_b2[oct][semi],
                                                           p->bp_a1[oct][semi], p->bp_a2[oct][semi],
                                                           bp_x1[oct][semi], bp_x2[oct][semi], 
                                                           bp_y1[oct][semi], bp_y2[oct][semi]);
                }
                
                // Apply saturation
                float sat_output = bypass_saturation ? combined_output : 
                                  tanhf(combined_output * sat_amount) * inv_internal_gain;
                
                // Apply LPF and mix
                float lpf_output = process_adaptive_lpf(sat_output, lpf_state, lpf2_state, lpf3_state, 
                                                       p->lpf_coeff, p->lpf_cutoff_smooth.current);
                
                buffer[frame] = input_sample * dry_gain + lpf_output * wet_gain;
                
                // Safety clipper
                if (fabsf(buffer[frame]) > 0.95f) {
                    buffer[frame] = tanhf(buffer[frame] * 0.7f) / 0.7f;
                }
            }
            break;
        }
        
        case 3: {
            // Three filters
            for (uint32_t frame = 0; frame < frames; ++frame) {
                float input_sample = buffer[frame];
                float combined_output = 0.0f;
                
                // Unroll the three filter operations
                for (int i = 0; i < 3; i++) {
                    int oct = p->active_filters[i].octave_idx;
                    int semi = p->active_filters[i].semitone_idx;
                    
                    combined_output += process_single_biquad(input_sample,
                                                           p->bp_b0[oct][semi], p->bp_b1[oct][semi], p->bp_b2[oct][semi],
                                                           p->bp_a1[oct][semi], p->bp_a2[oct][semi],
                                                           bp_x1[oct][semi], bp_x2[oct][semi], 
                                                           bp_y1[oct][semi], bp_y2[oct][semi]);
                }
                
                // Apply saturation
                float sat_output = bypass_saturation ? combined_output : 
                                  tanhf(combined_output * sat_amount) * inv_internal_gain;
                
                // Apply LPF and mix
                float lpf_output = process_adaptive_lpf(sat_output, lpf_state, lpf2_state, lpf3_state, 
                                                       p->lpf_coeff, p->lpf_cutoff_smooth.current);
                
                buffer[frame] = input_sample * dry_gain + lpf_output * wet_gain;
                
                // Safety clipper
                if (fabsf(buffer[frame]) > 0.95f) {
                    buffer[frame] = tanhf(buffer[frame] * 0.7f) / 0.7f;
                }
            }
            break;
        }
        
        case 4: {
            // Four filters
            for (uint32_t frame = 0; frame < frames; ++frame) {
                float input_sample = buffer[frame];
                float combined_output = 0.0f;
                
                // Unroll the four filter operations
                for (int i = 0; i < 4; i++) {
                    int oct = p->active_filters[i].octave_idx;
                    int semi = p->active_filters[i].semitone_idx;
                    
                    combined_output += process_single_biquad(input_sample,
                                                           p->bp_b0[oct][semi], p->bp_b1[oct][semi], p->bp_b2[oct][semi],
                                                           p->bp_a1[oct][semi], p->bp_a2[oct][semi],
                                                           bp_x1[oct][semi], bp_x2[oct][semi], 
                                                           bp_y1[oct][semi], bp_y2[oct][semi]);
                }
                
                // Apply saturation
                float sat_output = bypass_saturation ? combined_output : 
                                  tanhf(combined_output * sat_amount) * inv_internal_gain;
                
                // Apply LPF and mix
                float lpf_output = process_adaptive_lpf(sat_output, lpf_state, lpf2_state, lpf3_state, 
                                                       p->lpf_coeff, p->lpf_cutoff_smooth.current);
                
                buffer[frame] = input_sample * dry_gain + lpf_output * wet_gain;
                
                // Safety clipper
                if (fabsf(buffer[frame]) > 0.95f) {
                    buffer[frame] = tanhf(buffer[frame] * 0.7f) / 0.7f;
                }
            }
            break;
        }
        
        default: {
            // General case for 5+ filters
            for (uint32_t frame = 0; frame < frames; ++frame) {
                float input_sample = buffer[frame];
                float combined_output = 0.0f;
                
                // Process through active filters
                for (int i = 0; i < p->num_active_filters; i++) {
                    int oct = p->active_filters[i].octave_idx;
                    int semi = p->active_filters[i].semitone_idx;
                    
                    combined_output += process_single_biquad(input_sample,
                                                           p->bp_b0[oct][semi], p->bp_b1[oct][semi], p->bp_b2[oct][semi],
                                                           p->bp_a1[oct][semi], p->bp_a2[oct][semi],
                                                           bp_x1[oct][semi], bp_x2[oct][semi], 
                                                           bp_y1[oct][semi], bp_y2[oct][semi]);
                }
                
                // Apply saturation
                float sat_output = bypass_saturation ? combined_output : 
                                  tanhf(combined_output * sat_amount) * inv_internal_gain;
                
                // Apply LPF and mix
                float lpf_output = process_adaptive_lpf(sat_output, lpf_state, lpf2_state, lpf3_state, 
                                                       p->lpf_coeff, p->lpf_cutoff_smooth.current);
                
                buffer[frame] = input_sample * dry_gain + lpf_output * wet_gain;
                
                // Safety clipper
                if (fabsf(buffer[frame]) > 0.95f) {
                    buffer[frame] = tanhf(buffer[frame] * 0.7f) / 0.7f;
                }
            }
            break;
        }
    }
}

// Fixed adaptive LPF processing - proper cascaded filtering
float process_adaptive_lpf(float combined_output, float &lpf_state, float &lpf2_state, float &lpf3_state, float lpf_coeff, float lpf_cutoff_smooth) {
    float cutoff_normalized = (lpf_cutoff_smooth - 4000.0f) / 6000.0f;  // 0 at 4khz, 1 at 10khz
    cutoff_normalized = fmaxf(0.0f, fminf(1.0f, cutoff_normalized));
    
    // First pole (always active)
    lpf_state += lpf_coeff * (combined_output - lpf_state);
    float output = lpf_state;
    
    // Second pole (fade in above 4khz)
    if (cutoff_normalized > 0.0f) {
        lpf2_state += lpf_coeff * (output - lpf2_state);
        output = output * (1.0f - cutoff_normalized) + lpf2_state * cutoff_normalized;
    }
    
    // Third pole (fade in above 7khz) 
    if (cutoff_normalized > 0.5f) {
        float third_blend = (cutoff_normalized - 0.5f) * 2.0f;
        lpf3_state += lpf_coeff * (output - lpf3_state);
        output = output * (1.0f - third_blend) + lpf3_state * third_blend;
    }
    
    return output;
}

// Simplified oversampling functions (removed SIMD for better predictability)
void upsample_2x(const float *input, float *output, uint32_t input_frames, aa_filter_t *aa_filter) {
    // Simple scalar implementation - more predictable than SIMD for small buffers
    for (uint32_t i = 0; i < input_frames; i++) {
        output[i * 2] = input[i];
        output[i * 2 + 1] = 0.0f;
        
        output[i * 2] = process_aa_filter(output[i * 2], aa_filter);
        output[i * 2 + 1] = process_aa_filter(output[i * 2 + 1], aa_filter);
    }
}

void downsample_2x(const float *input, float *output, uint32_t input_frames, aa_filter_t *aa_filter) {
    uint32_t output_frames = input_frames / 2;
    
    // Simple scalar implementation
    for (uint32_t i = 0; i < output_frames; i++) {
        float sample1 = process_aa_filter(input[i * 2], aa_filter);
        float sample2 = process_aa_filter(input[i * 2 + 1], aa_filter);
        
        output[i] = (sample1 + sample2) * 0.5f;  // average both samples
    }
}
