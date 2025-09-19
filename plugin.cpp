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
                        p->root_note = fmax(ROOT_NOTE_MIN, fmin(ROOT_NOTE_MAX, param_event->value));
                        p->coefficients_need_update = true;
                        break;
                    case PARAM_BANDWIDTH:
                        p->bandwidth = fmax(BANDWIDTH_MIN, fmin(BANDWIDTH_MAX, param_event->value));
                        params_changed = true;
                        break;
                    case PARAM_ROOT_GAIN:
                        p->root_gain = fmax(ROOT_GAIN_MIN, fmin(ROOT_GAIN_MAX, param_event->value));
                        params_changed = true;
                        break;
                    case PARAM_LO_OCTAVES:
                        p->lo_octaves = fmax(LO_OCTAVES_MIN, fmin(LO_OCTAVES_MAX, param_event->value));
                        params_changed = true;
                        break;
                    case PARAM_HI_OCTAVES:
                        p->hi_octaves = fmax(HI_OCTAVES_MIN, fmin(HI_OCTAVES_MAX, param_event->value));
                        params_changed = true;
                        break;
                    case PARAM_HARMONICS:
                        p->harmonics = fmax(HARMONICS_MIN, fmin(HARMONICS_MAX, param_event->value));
                        params_changed = true;
                        break;
                    case PARAM_SATURATION:
                        p->saturation = fmax(SATURATION_MIN, fmin(SATURATION_MAX, param_event->value));
                        params_changed = true;
                        break;
                    case PARAM_LPF_CUTOFF:
                        p->lpf_cutoff = fmax(LPF_CUTOFF_MIN, fmin(LPF_CUTOFF_MAX, param_event->value));
                        params_changed = true;
                        break;
                    case PARAM_WET_BOOST:
                        p->wet_boost = fmax(WET_BOOST_MIN, fmin(WET_BOOST_MAX, param_event->value));
                        params_changed = true;
                        break;
                    case PARAM_DRY_WET:
                        p->dry_wet = fmax(DRY_WET_MIN, fmin(DRY_WET_MAX, param_event->value));
                        params_changed = true;
                        break;
                    case PARAM_SPREAD:
                        p->spread = fmax(SPREAD_MIN, fmin(SPREAD_MAX, param_event->value));
                        params_changed = true;
                        break;
                }
            }
        }
        
        // Only trigger smoothing if parameters actually changed
        if (params_changed) {
            trigger_parameter_smoothing(p);
        }
    }
    
    // check for valid inputs/outputs
    if (!process || process->audio_inputs_count == 0 || process->audio_outputs_count == 0) {
        return CLAP_PROCESS_CONTINUE;
    }
    
    const uint32_t nframes = process->frames_count;
    if (nframes == 0) return CLAP_PROCESS_CONTINUE;
    
    const clap_audio_buffer_t *input = &process->audio_inputs[0];
    const clap_audio_buffer_t *output = &process->audio_outputs[0];
    
    if (!input->data32 || !output->data32) {
        return CLAP_PROCESS_CONTINUE;
    }
    
    // Check if oversampling buffers are large enough
    uint32_t required_size = nframes * OVERSAMPLE_FACTOR;
    if (required_size > p->oversample_buffer_size || !p->oversample_buffer_L || !p->oversample_buffer_R) {
        return CLAP_PROCESS_CONTINUE;
    }
    
    // Process parameter smoothing once for the entire block
    process_parameter_smoothing(p, nframes);
    
    // check if we need to update filter coefficients (trust the flag)
    if (p->coefficients_need_update) {
        update_filter_coefficients(p);
    }
    
    // process audio with oversampling
    const uint32_t channels = input->channel_count < output->channel_count ? 
                             input->channel_count : output->channel_count;
    
    const uint32_t oversampled_frames = nframes * OVERSAMPLE_FACTOR;
    
    // Process left channel
    if (channels >= 1 && input->data32[0] && output->data32[0]) {
        // 1. Upsample input with proper interpolation
        upsample_2x(input->data32[0], p->oversample_buffer_L, nframes, &p->upsample_aa_L);
        
        // 2. Process oversampled signal
        process_channel_filters(p, p->oversample_buffer_L, oversampled_frames, true);
        
        // 3. Downsample and output
        downsample_2x(p->oversample_buffer_L, output->data32[0], oversampled_frames, &p->downsample_aa_L);
    }
    
    // Process right channel
    if (channels >= 2 && input->data32[1] && output->data32[1]) {
        // 1. Upsample input with proper interpolation
        upsample_2x(input->data32[1], p->oversample_buffer_R, nframes, &p->upsample_aa_R);
        
        // 2. Process oversampled signal
        process_channel_filters(p, p->oversample_buffer_R, oversampled_frames, false);
        
        // 3. Downsample and output
        downsample_2x(p->oversample_buffer_R, output->data32[1], oversampled_frames, &p->downsample_aa_R);
    }
    
    return CLAP_PROCESS_CONTINUE;
}

const void *minimal_get_extension(const clap_plugin *plugin, const char *id) {
    if (strcmp(id, CLAP_EXT_PARAMS) == 0) {
        return &plugin_params;
    }
    if (strcmp(id, CLAP_EXT_STATE) == 0) {
        return &plugin_state;
    }
    if (strcmp(id, CLAP_EXT_AUDIO_PORTS) == 0) {
        return &plugin_audio_ports;
    }
    return nullptr;
}

void minimal_on_main_thread(const clap_plugin *plugin) {
}

// plugin creation
static const clap_plugin_t *create_plugin(const clap_plugin_factory_t *factory,
                                          const clap_host_t *host,
                                          const char *plugin_id) {
    if (!plugin_id || strcmp(plugin_id, PLUGIN_ID) != 0) {
        return nullptr;
    }
    
    struct plugin_instance {
        clap_plugin_t plugin;
        minimal_plugin_t data;
    };
    
    plugin_instance *instance = (plugin_instance*)calloc(1, sizeof(plugin_instance));
    if (!instance) return nullptr;
    
    // setup plugin data
    instance->data.host = (clap_host_t*)host;
    instance->data.fine_tune = FINE_TUNE_DEFAULT;
    instance->data.root_note = ROOT_NOTE_DEFAULT;
    instance->data.bandwidth = BANDWIDTH_DEFAULT;
    instance->data.root_gain = ROOT_GAIN_DEFAULT;
    instance->data.lo_octaves = LO_OCTAVES_DEFAULT;
    instance->data.hi_octaves = HI_OCTAVES_DEFAULT;
    instance->data.harmonics = HARMONICS_DEFAULT;
    instance->data.saturation = SATURATION_DEFAULT;
    instance->data.lpf_cutoff = LPF_CUTOFF_DEFAULT;
    instance->data.wet_boost = WET_BOOST_DEFAULT;
    instance->data.dry_wet = DRY_WET_DEFAULT;
    instance->data.spread = SPREAD_DEFAULT;
    instance->data.sample_rate = 44100.0;
    instance->data.oversampled_rate = 44100.0 * OVERSAMPLE_FACTOR;
    instance->data.max_block_size = 0;  // will be set in activate()
    instance->data.lpf_state_L = 0.0f;
    instance->data.lpf_state_R = 0.0f;
    instance->data.lpf2_state_L = 0.0f;
    instance->data.lpf2_state_R = 0.0f;
    instance->data.lpf3_state_L = 0.0f;
    instance->data.lpf3_state_R = 0.0f;
    instance->data.lpf_coeff = 0.5f;
    
    // initialize oversampling buffers as null (will be allocated in activate)
    instance->data.oversample_buffer_L = nullptr;
    instance->data.oversample_buffer_R = nullptr;
    instance->data.oversample_buffer_size = 0;
    
    reset_all_filter_states(&instance->data);
    
    instance->data.last_fine_tune = -999.0;
    instance->data.last_root_note = -1.0;
    instance->data.last_bandwidth = -1.0;
    instance->data.last_lo_octaves = -1.0;
    instance->data.last_hi_octaves = -1.0;
    instance->data.last_harmonics = -1.0;
    instance->data.last_lpf_cutoff = -1.0;
    instance->data.last_root_gain = -1.0;  // ADDED
    instance->data.coefficients_need_update = true;
    
    // Initialize parameter smoothing
    update_parameter_smoothing(&instance->data);
    
    // setup plugin interface
    instance->plugin.desc = &plugin_desc;
    instance->plugin.plugin_data = &instance->data;
    instance->plugin.init = minimal_init;
    instance->plugin.destroy = minimal_destroy;
    instance->plugin.activate = minimal_activate;
    instance->plugin.deactivate = minimal_deactivate;
    instance->plugin.start_processing = minimal_start_processing;
    instance->plugin.stop_processing = minimal_stop_processing;
    instance->plugin.reset = minimal_reset;
    instance->plugin.process = minimal_process;
    instance->plugin.get_extension = minimal_get_extension;
    instance->plugin.on_main_thread = minimal_on_main_thread;
    
    return &instance->plugin;
}

// plugin factory
static uint32_t factory_get_plugin_count(const clap_plugin_factory_t *factory) {
    return 1;
}

static const clap_plugin_descriptor_t *factory_get_plugin_descriptor(const clap_plugin_factory_t *factory,
                                                                     uint32_t index) {
    return index == 0 ? &plugin_desc : nullptr;
}

static const clap_plugin_factory_t plugin_factory = {
    factory_get_plugin_count,
    factory_get_plugin_descriptor,
    create_plugin
};

// entry point factory function
static const void *get_factory(const char *factory_id) {
    if (strcmp(factory_id, CLAP_PLUGIN_FACTORY_ID) == 0) {
        return &plugin_factory;
    }
    return nullptr;
}

// main entry point
CLAP_EXPORT const clap_plugin_entry_t clap_entry = {
    CLAP_VERSION,
    [](const char *) -> bool { return true; },  // init
    []() -> void {},                             // deinit  
    get_factory
};
