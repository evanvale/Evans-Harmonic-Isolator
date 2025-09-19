#include "plugin.h"
#include <cstring>
#include <cstdio>
#include <stdlib.h>
#include <math.h>

// plugin features
static const char *features[] = {
    CLAP_PLUGIN_FEATURE_AUDIO_EFFECT,
    nullptr
};

// plugin descriptor
static const clap_plugin_descriptor_t plugin_desc = {
    CLAP_VERSION,
    PLUGIN_ID,
    PLUGIN_NAME,
    PLUGIN_VENDOR,
    nullptr,        // url
    nullptr,        // manual_url  
    nullptr,        // support_url
    PLUGIN_VERSION,
    PLUGIN_DESC,
    features
};

// plugin lifecycle implementations
bool minimal_init(const clap_plugin *plugin) {
    minimal_plugin_t *p = (minimal_plugin_t*)plugin->plugin_data;
    
    // initialize parameter values
    p->fine_tune = FINE_TUNE_DEFAULT;
    p->root_note = ROOT_NOTE_DEFAULT;
    p->bandwidth = BANDWIDTH_DEFAULT;
    p->root_gain = ROOT_GAIN_DEFAULT;
    p->lo_octaves = LO_OCTAVES_DEFAULT;
    p->hi_octaves = HI_OCTAVES_DEFAULT;
    p->harmonics = HARMONICS_DEFAULT;
    p->saturation = SATURATION_DEFAULT;
    p->lpf_cutoff = LPF_CUTOFF_DEFAULT;
    p->wet_boost = WET_BOOST_DEFAULT;
    p->dry_wet = DRY_WET_DEFAULT;
    p->spread = SPREAD_DEFAULT;
    
    // initialize DSP state
    p->sample_rate = 44100.0;
    p->oversampled_rate = p->sample_rate * OVERSAMPLE_FACTOR;
    p->max_block_size = 0;  // will be set in activate()
    p->lpf_state_L = 0.0f;
    p->lpf_state_R = 0.0f;
    p->lpf2_state_L = 0.0f;
    p->lpf2_state_R = 0.0f;
    p->lpf3_state_L = 0.0f;
    p->lpf3_state_R = 0.0f;
    p->lpf_coeff = 0.5f;
    
    // initialize oversampling buffers
    p->oversample_buffer_L = nullptr;
    p->oversample_buffer_R = nullptr;
    p->oversample_buffer_size = 0;
    
    // initialize all filter states
    reset_all_filter_states(p);
    
    // initialize parameter tracking
    p->last_fine_tune = -999.0;  // Make sure it updates on first run
    p->last_root_note = -1.0;
    p->last_bandwidth = -1.0;
    p->last_lo_octaves = -1.0;
    p->last_hi_octaves = -1.0;
    p->last_harmonics = -1.0;
    p->last_lpf_cutoff = -1.0;
    p->last_root_gain = -1.0;  // ADDED
    p->coefficients_need_update = true;
    
    // initialize parameter smoothing
    update_parameter_smoothing(p);
    
    // Ensure all smoothed values start at their parameter values
    p->bandwidth_smooth.current = p->bandwidth;
    p->root_gain_smooth.current = p->root_gain;
    p->lo_octaves_smooth.current = p->lo_octaves;
    p->hi_octaves_smooth.current = p->hi_octaves;
    p->harmonics_smooth.current = p->harmonics;
    p->saturation_smooth.current = p->saturation;
    p->lpf_cutoff_smooth.current = p->lpf_cutoff;
    p->wet_boost_smooth.current = p->wet_boost;
    p->dry_wet_smooth.current = p->dry_wet;
    p->spread_smooth.current = p->spread;
    
    // calculate initial coefficients to avoid hiccup on first parameter change
    update_filter_coefficients(p);
    
    return true;
}

void minimal_destroy(const clap_plugin *plugin) {
    if (plugin) {
        minimal_plugin_t *p = (minimal_plugin_t*)plugin->plugin_data;
        cleanup_oversampling(p);
        free((void*)plugin);
    }
}

bool minimal_activate(const clap_plugin *plugin, double sample_rate, uint32_t min_frames, uint32_t max_frames) {
    minimal_plugin_t *p = (minimal_plugin_t*)plugin->plugin_data;
    
    p->sample_rate = sample_rate;
    p->oversampled_rate = sample_rate * OVERSAMPLE_FACTOR;
    p->max_block_size = max_frames;  // Store max_frames for proper buffer allocation
    
    // Initialize oversampling with proper buffer size
    init_oversampling(p);
    
    // Initialize parameter smoothing with new sample rate
    update_parameter_smoothing(p);
    
    reset_all_filter_states(p);
    p->coefficients_need_update = true;
    p->lpf_state_L = 0.0f;
    p->lpf_state_R = 0.0f;
    p->lpf2_state_L = 0.0f;
    p->lpf2_state_R = 0.0f;
    p->lpf3_state_L = 0.0f;
    p->lpf3_state_R = 0.0f;
    
    // sync tracking variables to avoid hiccup on first parameter change
    p->last_fine_tune = p->fine_tune;
    p->last_root_note = p->root_note;
    p->last_bandwidth = p->bandwidth_smooth.current;
    p->last_lo_octaves = p->lo_octaves_smooth.current;
    p->last_hi_octaves = p->hi_octaves_smooth.current;
    p->last_harmonics = p->harmonics_smooth.current;
    p->last_lpf_cutoff = p->lpf_cutoff_smooth.current;
    p->last_root_gain = p->root_gain_smooth.current;
    
    return true;
}

void minimal_deactivate(const clap_plugin *plugin) {
    minimal_plugin_t *p = (minimal_plugin_t*)plugin->plugin_data;
    cleanup_oversampling(p);
}

bool minimal_start_processing(const clap_plugin *plugin) {
    return true;
}

void minimal_stop_processing(const clap_plugin *plugin) {
}

void minimal_reset(const clap_plugin *plugin) {
    minimal_plugin_t *p = (minimal_plugin_t*)plugin->plugin_data;
    reset_all_filter_states(p);
    p->lpf_state_L = 0.0f;
    p->lpf_state_R = 0.0f;
    p->lpf2_state_L = 0.0f;
    p->lpf2_state_R = 0.0f;
    p->lpf3_state_L = 0.0f;
    p->lpf3_state_R = 0.0f;
}

clap_process_status minimal_process(const clap_plugin *plugin, const clap_process_t *process) {
    minimal_plugin_t *p = (minimal_plugin_t*)plugin->plugin_data;
    
    // Handle parameter events
    if (process->in_events) {
        uint32_t event_count = process->in_events->size(process->in_events);
        bool params_changed = false;
        
        for (uint32_t i = 0; i < event_count; ++i) {
            const clap_event_header_t *header = process->in_events->get(process->in_events, i);
            if (header->type == CLAP_EVENT_PARAM_VALUE) {
                const clap_event_param_value_t *param_event = (const clap_event_param_value_t*)header;
                
                switch (param_event->param_id) {
                    case PARAM_FINE_TUNE:
                        p->fine_tune = fmax(FINE_TUNE_MIN, fmin(FINE_TUNE_MAX, param_event->value));
                        p->coefficients_need_update = true;
                        break;
                    case PARAM_ROOT_NOTE:
                        p->root_note = fmax(ROOT
