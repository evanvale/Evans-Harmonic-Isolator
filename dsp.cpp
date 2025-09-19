#include "plugin.h"
#include <cmath>
#include <cstring>
#include <cstdlib>

// Cross-platform SIMD support
#ifdef __APPLE__
    // Apple universal binary support
    #if defined(__x86_64__) || defined(__i386__)
        #include <immintrin.h>  // x86 SIMD
    #elif defined(__aarch64__) || defined(__arm64__)
        // ARM64: disable x86 SIMD for now
        #undef __SSE__
    #endif
#else
    // Non-Apple platforms - assume x86
    #include <immintrin.h>
#endif

static const double PI = 3.14159265358979323846;

// Tunable constants - UPDATED FOR OVERSAMPLING
static const double MIN_FREQ = 50.0;              
static const double MAX_FREQ = 5000.0;
static const double UPPER_NYQUIST_LIMIT = 0.45;   // CHANGED: was 0.475
static const double MAX_Q = 50.0;                 
static const double MIN_Q = 2.0;                  
static const double LOW_FREQ_Q_KNEE = 140.0;      
static const double ROOT_Q_SCALE_DIV = 300.0;     
static const double OCTAVE_TAPER_BASE = 0.8;      
static const double MIN_TAPER = 0.035;            

static const double GAIN_SOFT_Q_REF = 25.0;       
static const double GAIN_SOFT_STRENGTH = 0.03;    

// Musical harmonic intervals
static const int HARMONIC_5TH_SEMITONES = 7;      
static const int HARMONIC_9TH_SEMITONES = 14;     
static const int HARMONIC_13TH_SEMITONES = 21;    

static const double HARMONIC_5TH_GAIN = 1.0;      
static const double HARMONIC_9TH_GAIN = 0.7;      
static const double HARMONIC_13TH_GAIN = 0.4;     

// Enable flush-to-zero for denormal protection
static void enable_flush_to_zero() {
    #if defined(__SSE__)
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    #endif
}

// Simplified anti-aliasing filter coefficients for 2x oversampling
// 2nd order Butterworth @ 0.45 Nyquist (much simpler than 8th order elliptic)
static const float AA_FILTER_COEFFS[2][5] = {
    {0.0976f, 0.1953f, 0.0976f, -1.0177f, 0.4082f},  // 2nd order butterworth
    {0.0976f, 0.1953f, 0.0976f, -1.0177f, 0.4082f}   // 2nd order butterworth
};

static inline double smooth_gain_comp(double Q) {
    double over = fmax(0.0, Q - GAIN_SOFT_Q_REF);
    return 1.0 / (1.0 + GAIN_SOFT_STRENGTH * over * over);
}

// Smooth bass boost function to eliminate discontinuities
static inline double smooth_bass_boost(double freq) {
    if (freq >= 250.0) return 1.0;
    
    // Use smooth curves instead of steps
    if (freq <= 60.0) {
        // Very low frequencies: smooth curve from 4.0 to current
        double t = freq / 60.0;
        t = t * t; // square for smoother curve
        return 4.0 * (0.5 + 0.5 * t) + 4.0 * (1.0 - t);  // 4.0 to 8.0 range
    } else if (freq <= 100.0) {
        // Sub-bass range: 6.0 to 4.0
        double t = (freq - 60.0) / 40.0;
        t = 0.5 * (1.0 - cos(PI * t)); // smooth cosine interpolation
        return 6.0 * (1.0 - t) + 4.0 * t;
    } else if (freq <= 150.0) {
        // Low bass range: 4.0 to 2.5
        double t = (freq - 100.0) / 50.0;
        t = 0.5 * (1.0 - cos(PI * t));
        return 4.0 * (1.0 - t) + 2.5 * t;
    } else {
        // Upper bass range: 2.5 to 1.0
        double t = (freq - 150.0) / 100.0;
        t = 0.5 * (1.0 - cos(PI * t));
        return 2.5 * (1.0 - t) + 1.0 * t;
    }
}

void init_oversampling(minimal_plugin_t *p) {
    // Enable flush-to-zero for denormal protection
    enable_flush_to_zero();
    
    // Use max_block_size from activate() for proper buffer allocation
    p->oversample_buffer_size = p->max_block_size * OVERSAMPLE_FACTOR;
    
    // Align buffers for SIMD operations
    #if defined(__SSE__)
    p->oversample_buffer_L = (float*)_mm_malloc(p->oversample_buffer_size * sizeof(float), 16);
    p->oversample_buffer_R = (float*)_mm_malloc(p->oversample_buffer_size * sizeof(float), 16);
    #else
    p->oversample_buffer_L = (float*)malloc(p->oversample_buffer_size * sizeof(float));
    p->oversample_buffer_R = (float*)malloc(p->oversample_buffer_size * sizeof(float));
    #endif
    
    if (!p->oversample_buffer_L || !p->oversample_buffer_R) {
        cleanup_oversampling(p);
        return;  // allocation failed
    }
    
    // Clear buffers
    memset(p->oversample_buffer_L, 0, p->oversample_buffer_size * sizeof(float));
    memset(p->oversample_buffer_R, 0, p->oversample_buffer_size * sizeof(float));
    
    // Initialize simplified AA filter coefficients (2 biquads instead of 4)
    for (int i = 0; i < 2; i++) {
        // Upsample filters
        p->upsample_aa_L.b0[i] = AA_FILTER_COEFFS[i][0];
        p->upsample_aa_L.b1[i] = AA_FILTER_COEFFS[i][1];
        p->upsample_aa_L.b2[i] = AA_FILTER_COEFFS[i][2];
        p->upsample_aa_L.a1[i] = AA_FILTER_COEFFS[i][3];
        p->upsample_aa_L.a2[i] = AA_FILTER_COEFFS[i][4];
        
        p->upsample_aa_R.b0[i] = AA_FILTER_COEFFS[i][0];
        p->upsample_aa_R.b1[i] = AA_FILTER_COEFFS[i][1];
        p->upsample_aa_R.b2[i] = AA_FILTER_COEFFS[i][2];
        p->upsample_aa_R.a1[i] = AA_FILTER_COEFFS[i][3];
        p->upsample_aa_R.a2[i] = AA_FILTER_COEFFS[i][4];
        
        // Downsample filters (same coefficients)
        p->downsample_aa_L.b0[i] = AA_FILTER_COEFFS[i][0];
        p->downsample_aa_L.b1[i] = AA_FILTER_COEFFS[i][1];
        p->downsample_aa_L.b2[i] = AA_FILTER_COEFFS[i][2];
        p->downsample_aa_L.a1[i] = AA_FILTER_COEFFS[i][3];
        p->downsample_aa_L.a2[i] = AA_FILTER_COEFFS[i][4];
        
        p->downsample_aa_R.b0[i] = AA_FILTER_COEFFS[i][0];
        p->downsample_aa_R.b1[i] = AA_FILTER_COEFFS[i][1];
        p->downsample_aa_R.b2[i] = AA_FILTER_COEFFS[i][2];
        p->downsample_aa_R.a1[i] = AA_FILTER_COEFFS[i][3];
        p->downsample_aa_R.a2[i] = AA_FILTER_COEFFS[i][4];
    }
    
    // Zero out unused filter stages
    for (int i = 2; i < 4; i++) {
        p->upsample_aa_L.b0[i] = p->upsample_aa_L.b1[i] = p->upsample_aa_L.b2[i] = 0.0f;
        p->upsample_aa_L.a1[i] = p->upsample_aa_L.a2[i] = 0.0f;
        p->upsample_aa_R.b0[i] = p->upsample_aa_R.b1[i] = p->upsample_aa_R.b2[i] = 0.0f;
        p->upsample_aa_R.a1[i] = p->upsample_aa_R.a2[i] = 0.0f;
        p->downsample_aa_L.b0[i] = p->downsample_aa_L.b1[i] = p->downsample_aa_L.b2[i] = 0.0f;
        p->downsample_aa_L.a1[i] = p->downsample_aa_L.a2[i] = 0.0f;
        p->downsample_aa_R.b0[i] = p->downsample_aa_R.b1[i] = p->downsample_aa_R.b2[i] = 0.0f;
        p->downsample_aa_R.a1[i] = p->downsample_aa_R.a2[i] = 0.0f;
    }
    
    reset_aa_filters(p);
}

void cleanup_oversampling(minimal_plugin_t *p) {
    if (p->oversample_buffer_L) {
        #if defined(__SSE__)
        _mm_free(p->oversample_buffer_L);
        #else
        free(p->oversample_buffer_L);
        #endif
        p->oversample_buffer_L = nullptr;
    }
    if (p->oversample_buffer_R) {
        #if defined(__SSE__)
        _mm_free(p->oversample_buffer_R);
        #else
        free(p->oversample_buffer_R);
        #endif
        p->oversample_buffer_R = nullptr;
    }
    p->oversample_buffer_size = 0;
}

void reset_aa_filters(minimal_plugin_t *p) {
    for (int i = 0; i < 4; i++) {
        p->upsample_aa_L.x1[i] = p->upsample_aa_L.x2[i] = 0.0f;
        p->upsample_aa_L.y1[i] = p->upsample_aa_L.y2[i] = 0.0f;
        p->upsample_aa_R.x1[i] = p->upsample_aa_R.x2[i] = 0.0f;
        p->upsample_aa_R.y1[i] = p->upsample_aa_R.y2[i] = 0.0f;
        
        p->downsample_aa_L.x1[i] = p->downsample_aa_L.x2[i] = 0.0f;
        p->downsample_aa_L.y1[i] = p->downsample_aa_L.y2[i] = 0.0f;
        p->downsample_aa_R.x1[i] = p->downsample_aa_R.x2[i] = 0.0f;
        p->downsample_aa_R.y1[i] = p->downsample_aa_R.y2[i] = 0.0f;
    }
}

float process_aa_filter(float input, aa_filter_t *filter) {
    float output = input;
    
    // Process through only 2 biquads instead of 4 (4th order instead of 8th)
    for (int i = 0; i < 2; i++) {
        float temp = filter->b0[i] * output + filter->b1[i] * filter->x1[i] + filter->b2[i] * filter->x2[i]
                    - filter->a1[i] * filter->y1[i] - filter->a2[i] * filter->y2[i];
        
        filter->x2[i] = filter->x1[i];
        filter->x1[i] = output;
        filter->y2[i] = filter->y1[i];
        filter->y1[i] = temp;
        
        output = temp;
    }
    
    return output;
}

void reset_all_filter_states(minimal_plugin_t *p) {
    for (int oct = 0; oct < NUM_OCTAVES; oct++) {
        for (int semi = 0; semi < NUM_SEMITONES; semi++) {
            p->bp_x1_L[oct][semi] = p->bp_x2_L[oct][semi] = 0.0f;
            p->bp_y1_L[oct][semi] = p->bp_y2_L[oct][semi] = 0.0f;
            p->bp_x1_R[oct][semi] = p->bp_x2_R[oct][semi] = 0.0f;
            p->bp_y1_R[oct][semi] = p->bp_y2_R[oct][semi] = 0.0f;
        }
    }
    // Also reset LPF states
    p->lpf_state_L = 0.0f;
    p->lpf_state_R = 0.0f;
    p->lpf2_state_L = 0.0f;
    p->lpf2_state_R = 0.0f;
    p->lpf3_state_L = 0.0f;
    p->lpf3_state_R = 0.0f;
    
    // Reset active filter tracking
    p->num_active_filters = 0;
    
    // Reset AA filters
    reset_aa_filters(p);
}

// Optimized coefficient calculation using direct trig calls (no lookup tables)
static inline void set_bp_coeffs_constant_skirt_optimized(
    float &b0, float &b1, float &b2,
    float &a1, float &a2,
    double omega, double Q, double gain)
{
    // Use direct sinf/cosf - modern CPUs make these very fast
    const float sin_omega = sinf(omega);
    const float cos_omega = cosf(omega);
    const double alpha = sin_omega / (2.0 * Q);

    const double a0 = 1.0 + alpha;
    const double inv_a0 = 1.0 / a0;

    double B0 = sin_omega / (2.0 * Q);
    double B1 = 0.0;
    double B2 = -sin_omega / (2.0 * Q);
    double A1 = -2.0 * cos_omega;
    double A2 = 1.0 - alpha;

    b0 = static_cast<float>( (B0 * gain) * inv_a0 );
    b1 = static_cast<float>( (B1 * gain) * inv_a0 );
    b2 = static_cast<float>( (B2 * gain) * inv_a0 );
    a1 = static_cast<float>( A1 * inv_a0 );
    a2 = static_cast<float>( A2 * inv_a0 );
}

void update_parameter_smoothing(minimal_plugin_t *p) {
    // Calculate smoothing coefficient based on sample rate
    // 5ms smoothing time at base sample rate (not oversampled!)
    float smooth_time_samples = (5.0f / 1000.0f) * p->sample_rate;
    p->smooth_coeff = 1.0f / smooth_time_samples;
    
    // Initialize smoothed values to current parameter values
    p->bandwidth_smooth.current = p->bandwidth;
    p->bandwidth_smooth.target = p->bandwidth;
    p->bandwidth_smooth.active = false;
    
    p->root_gain_smooth.current = p->root_gain;
    p->root_gain_smooth.target = p->root_gain;
    p->root_gain_smooth.active = false;
    
    p->lo_octaves_smooth.current = p->lo_octaves;
    p->lo_octaves_smooth.target = p->lo_octaves;
    p->lo_octaves_smooth.active = false;
    
    p->hi_octaves_smooth.current = p->hi_octaves;
    p->hi_octaves_smooth.target = p->hi_octaves;
    p->hi_octaves_smooth.active = false;
    
    p->harmonics_smooth.current = p->harmonics;
    p->harmonics_smooth.target = p->harmonics;
    p->harmonics_smooth.active = false;
    
    p->saturation_smooth.current = p->saturation;
    p->saturation_smooth.target = p->saturation;
    p->saturation_smooth.active = false;
    
    p->lpf_cutoff_smooth.current = p->lpf_cutoff;
    p->lpf_cutoff_smooth.target = p->lpf_cutoff;
    p->lpf_cutoff_smooth.active = false;
    
    p->wet_boost_smooth.current = p->wet_boost;
    p->wet_boost_smooth.target = p->wet_boost;
    p->wet_boost_smooth.active = false;
    
    p->dry_wet_smooth.current = p->dry_wet;
    p->dry_wet_smooth.target = p->dry_wet;
    p->dry_wet_smooth.active = false;
    
    p->spread_smooth.current = p->spread;
    p->spread_smooth.target = p->spread;
    p->spread_smooth.active = false;
    
    p->any_smoothing_active = false;
}

// Trigger smoothing for a specific parameter
void trigger_parameter_smoothing(minimal_plugin_t *p) {
    const float epsilon = 0.0001f;
    
    // Check each parameter and activate smoothing if needed
    if (fabs(p->bandwidth - p->bandwidth_smooth.target) > epsilon) {
        p->bandwidth_smooth.target = p->bandwidth;
        p->bandwidth_smooth.active = true;
        p->any_smoothing_active = true;
    }
    
    if (fabs(p->root_gain - p->root_gain_smooth.target) > epsilon) {
        p->root_gain_smooth.target = p->root_gain;
        p->root_gain_smooth.active = true;
        p->any_smoothing_active = true;
    }
    
    if (fabs(p->lo_octaves - p->lo_octaves_smooth.target) > epsilon) {
        p->lo_octaves_smooth.target = p->lo_octaves;
        p->lo_octaves_smooth.active = true;
        p->any_smoothing_active = true;
    }
    
    if (fabs(p->hi_octaves - p->hi_octaves_smooth.target) > epsilon) {
        p->hi_octaves_smooth.target = p->hi_octaves;
        p->hi_octaves_smooth.active = true;
        p->any_smoothing_active = true;
    }
    
    if (fabs(p->harmonics - p->harmonics_smooth.target) > epsilon) {
        p->harmonics_smooth.target = p->harmonics;
        p->harmonics_smooth.active = true;
        p->any_smoothing_active = true;
    }
    
    if (fabs(p->saturation - p->saturation_smooth.target) > epsilon) {
        p->saturation_smooth.target = p->saturation;
        p->saturation_smooth.active = true;
        p->any_smoothing_active = true;
    }
    
    if (fabs(p->lpf_cutoff - p->lpf_cutoff_smooth.target) > epsilon) {
        p->lpf_cutoff_smooth.target = p->lpf_cutoff;
        p->lpf_cutoff_smooth.active = true;
        p->any_smoothing_active = true;
    }
    
    if (fabs(p->wet_boost - p->wet_boost_smooth.target) > epsilon) {
        p->wet_boost_smooth.target = p->wet_boost;
        p->wet_boost_smooth.active = true;
        p->any_smoothing_active = true;
    }
    
    if (fabs(p->dry_wet - p->dry_wet_smooth.target) > epsilon) {
        p->dry_wet_smooth.target = p->dry_wet;
        p->dry_wet_smooth.active = true;
        p->any_smoothing_active = true;
    }
    
    if (fabs(p->spread - p->spread_smooth.target) > epsilon) {
        p->spread_smooth.target = p->spread;
        p->spread_smooth.active = true;
        p->any_smoothing_active = true;
    }
}

// Process active smoothing - call once per process block
void process_parameter_smoothing(minimal_plugin_t *p, uint32_t frames) {
    // Early exit if nothing needs smoothing
    if (!p->any_smoothing_active) return;
    
    const float smooth_threshold = 0.001f;
    bool still_active = false;
    
    // Helper lambda to smooth a single parameter (batched every 32 samples)
    auto smooth_param_batched = [&](param_smooth_t &param) {
        if (!param.active) return;
        
        // Update in 32-sample chunks for better performance
        uint32_t batch_size = 32;
        uint32_t full_batches = frames / batch_size;
        uint32_t remainder = frames % batch_size;
        
        // Process full batches
        for (uint32_t batch = 0; batch < full_batches; batch++) {
            for (uint32_t i = 0; i < batch_size; i++) {
                param.current += (param.target - param.current) * p->smooth_coeff;
            }
        }
        
        // Process remainder samples
        for (uint32_t i = 0; i < remainder; i++) {
            param.current += (param.target - param.current) * p->smooth_coeff;
        }
        
        // Check if we're close enough to deactivate
        if (fabsf(param.current - param.target) < smooth_threshold) {
            param.current = param.target;
            param.active = false;
        } else {
            still_active = true;
        }
    };
    
    // Helper lambda for per-sample smoothing (only for critical parameters)
    auto smooth_param_per_sample = [&](param_smooth_t &param) {
        if (!param.active) return;
        
        // Apply exponential smoothing for the number of frames
        for (uint32_t i = 0; i < frames; i++) {
            param.current += (param.target - param.current) * p->smooth_coeff;
        }
        
        // Check if we're close enough to deactivate
        if (fabsf(param.current - param.target) < smooth_threshold) {
            param.current = param.target;
            param.active = false;
        } else {
            still_active = true;
        }
    };
    
    // Keep 32-sample smoothing for dry_wet to avoid zipper noise
    auto smooth_param_32 = [&](param_smooth_t &param) {
        if (!param.active) return;
        
        // Update in 32-sample chunks
        uint32_t batch_size = 32;
        uint32_t full_batches = frames / batch_size;
        uint32_t remainder = frames % batch_size;
        
        // Process full batches
        for (uint32_t batch = 0; batch < full_batches; batch++) {
            for (uint32_t i = 0; i < batch_size; i++) {
                param.current += (param.target - param.current) * p->smooth_coeff;
            }
        }
        
        // Process remainder samples
        for (uint32_t i = 0; i < remainder; i++) {
            param.current += (param.target - param.current) * p->smooth_coeff;
        }
        
        // Check if we're close enough to deactivate
        if (fabsf(param.current - param.target) < smooth_threshold) {
            param.current = param.target;
            param.active = false;
        } else {
            still_active = true;
        }
    };
    
    // Process parameters with 64-sample batched smoothing (less critical)
    smooth_param_batched(p->bandwidth_smooth);
    smooth_param_batched(p->root_gain_smooth);
    smooth_param_batched(p->lo_octaves_smooth);
    smooth_param_batched(p->hi_octaves_smooth);
    smooth_param_batched(p->harmonics_smooth);
    smooth_param_batched(p->saturation_smooth);
    smooth_param_batched(p->lpf_cutoff_smooth);
    smooth_param_batched(p->wet_boost_smooth);
    smooth_param_batched(p->spread_smooth);
    
    // Use 32-sample smoothing for dry_wet to avoid zipper noise
    smooth_param_32(p->dry_wet_smooth);
    
    // Update global active flag
    p->any_smoothing_active = still_active;
    
    // Check if smoothed values affecting coefficients have changed enough
    if (fabs(p->bandwidth_smooth.current - p->last_bandwidth) > smooth_threshold ||
        fabs(p->lo_octaves_smooth.current - p->last_lo_octaves) > smooth_threshold ||
        fabs(p->hi_octaves_smooth.current - p->last_hi_octaves) > smooth_threshold ||
        fabs(p->harmonics_smooth.current - p->last_harmonics) > smooth_threshold ||
        fabs(p->lpf_cutoff_smooth.current - p->last_lpf_cutoff) > smooth_threshold ||
        fabs(p->root_gain_smooth.current - p->last_root_gain) > smooth_threshold) {
        p->coefficients_need_update = true;
    }
}

void update_filter_coefficients(minimal_plugin_t *p) {
    // Apply fine tuning to root frequency
    double target_freq = semitone_to_freq(p->root_note + p->fine_tune / 100.0);

    double Q = 2.0 + (1.0 - p->bandwidth_smooth.current) * (MAX_Q - 2.0);
    if (Q < MIN_Q) Q = MIN_Q;

    const double base_gain_comp = smooth_gain_comp(Q);

    // Clear all filter coefficients and reset active filter tracking
    for (int oct = 0; oct < NUM_OCTAVES; oct++) {
        for (int semi = 0; semi < NUM_SEMITONES; semi++) {
            p->bp_b0[oct][semi] = 0.0f;
            p->bp_b1[oct][semi] = 0.0f;
            p->bp_b2[oct][semi] = 0.0f;
            p->bp_a1[oct][semi] = 0.0f;
            p->bp_a2[oct][semi] = 0.0f;
        }
    }
    
    p->num_active_filters = 0;

    // Cache common calculations
    const double inv_oversampled_rate = 1.0 / p->oversampled_rate;
    const double two_pi = 2.0 * PI;
    
    // Lambda for applying a single filter with active filter tracking
    auto apply_filter = [&](double freq, double q, double gain) {
        if (p->num_active_filters >= 16) return; // safety limit
        if (freq < MIN_FREQ || freq >= p->oversampled_rate * UPPER_NYQUIST_LIMIT || freq > MAX_FREQ) return;

        const double omega = two_pi * freq * inv_oversampled_rate;
        
        // Use the next available filter slot
        int oct_idx = p->num_active_filters;
        int semi_idx = 0;
        
        set_bp_coeffs_constant_skirt_optimized(
            p->bp_b0[oct_idx][semi_idx], p->bp_b1[oct_idx][semi_idx], p->bp_b2[oct_idx][semi_idx],
            p->bp_a1[oct_idx][semi_idx], p->bp_a2[oct_idx][semi_idx],
            omega, q, gain
        );
        
        // Track this as an active filter
        p->active_filters[p->num_active_filters].octave_idx = oct_idx;
        p->active_filters[p->num_active_filters].semitone_idx = semi_idx;
        p->num_active_filters++;
    };

    // Main root filter
    if (target_freq > MIN_FREQ && target_freq < p->oversampled_rate * 0.45) {
        double main_q_factor = fmin(1.0, fmax(0.6, target_freq / ROOT_Q_SCALE_DIV));
        double main_q = fmax(MIN_Q, Q * main_q_factor);

        if (target_freq < LOW_FREQ_Q_KNEE) {
            double t = target_freq / LOW_FREQ_Q_KNEE;
            double ease = 0.55 + 0.45 * t;
            main_q = fmax(MIN_Q, main_q * ease);
        }

        double root_gain = base_gain_comp * p->root_gain_smooth.current * smooth_bass_boost(target_freq);
        apply_filter(target_freq, main_q, root_gain);
    }

    // Musical harmonics (5th, 9th, 13th)
    if (p->harmonics_smooth.current > 0.0) {
        struct Harmonic { int semitones; double gain_mult; double q_mult; };
        Harmonic harmonics[] = {
            { HARMONIC_5TH_SEMITONES, HARMONIC_5TH_GAIN, 1.1 },
            { HARMONIC_9TH_SEMITONES, HARMONIC_9TH_GAIN, 1.15 },
            { HARMONIC_13TH_SEMITONES, HARMONIC_13TH_GAIN, 1.2 }
        };

        for (const auto& h : harmonics) {
            if (p->num_active_filters >= 16) break;
            
            double freq = target_freq * pow(2.0, h.semitones / 12.0);
            double harm_q = Q * h.q_mult;

            double gain = base_gain_comp * p->harmonics_smooth.current * h.gain_mult;
            apply_filter(freq, harm_q, gain);
        }
    }

    // High octaves (above root)
    if (p->hi_octaves_smooth.current > 0.0) {
        for (int oct = 1; oct <= 4 && p->num_active_filters < 16; oct++) {
            const double oct_freq = target_freq * pow(2.0, oct);
            double taper = p->hi_octaves_smooth.current * pow(OCTAVE_TAPER_BASE, oct);
            if (taper < MIN_TAPER) continue;

            double gain = base_gain_comp * taper;
            apply_filter(oct_freq, Q, gain);
        }
    }

// Low octaves (below root) - with enhanced bass compensation
    if (p->lo_octaves_smooth.current > 0.0) {
        for (int oct = 1; oct <= 3 && p->num_active_filters < NUM_OCTAVES; oct++) {
            const double oct_freq = target_freq / pow(2.0, oct);
            if (oct_freq < MIN_FREQ) continue;
            
            double taper = p->lo_octaves_smooth.current * pow(OCTAVE_TAPER_BASE, oct);
            if (taper < MIN_TAPER) continue;

            double q = Q;
            if (oct_freq < LOW_FREQ_Q_KNEE) {
                double t = oct_freq / LOW_FREQ_Q_KNEE;
                q = fmax(MIN_Q, Q * (0.65 + 0.35 * t));
            }

            // Use smooth bass boost function
            double gain = base_gain_comp * taper * smooth_bass_boost(oct_freq);
            
            // Also widen Q significantly for fatter bass
            if (oct_freq < 250.0) {
                q *= 0.6;
                if (q < MIN_Q) q = MIN_Q;
            }

            apply_filter(oct_freq, q, gain);
        }
    }

    // One-pole LPF coefficient (now at oversampled rate)
    p->lpf_coeff = 1.0f - expf(-2.0f * M_PI * p->lpf_cutoff_smooth.current / p->oversampled_rate);

    // Tracking variables - track smoothed values, not raw ones
    p->last_fine_tune = p->fine_tune;  // This one is not smoothed
    p->last_root_note = p->root_note;  // This one is not smoothed  
    p->last_bandwidth = p->bandwidth_smooth.current;
    p->last_lo_octaves = p->lo_octaves_smooth.current;
    p->last_hi_octaves = p->hi_octaves_smooth.current;
    p->last_harmonics = p->harmonics_smooth.current;
    p->last_lpf_cutoff = p->lpf_cutoff_smooth.current;
    p->last_root_gain = p->root_gain_smooth.current;
    p->coefficients_need_update = false;
}
