#include "plugin.h"
#include <cstring>
#include <cstdio>
#include <cmath>

static uint32_t params_count(const clap_plugin *plugin) {
    return PARAM_COUNT;  // Should be 12 now!
}

static bool params_get_info(const clap_plugin *plugin, uint32_t param_index, clap_param_info_t *param_info) {
    if (!param_info) return false;
    // zero everything first to avoid junk memory
    clap_param_info_t info{};
    
    switch (param_index) {
        case PARAM_FINE_TUNE:
            info.id = PARAM_FINE_TUNE;
            info.flags = CLAP_PARAM_IS_AUTOMATABLE;
            info.min_value = FINE_TUNE_MIN;
            info.max_value = FINE_TUNE_MAX;
            info.default_value = FINE_TUNE_DEFAULT;
            strncpy(info.name, "Fine Tune", CLAP_NAME_SIZE-1);
            strncpy(info.module, "", CLAP_PATH_SIZE-1);
            break;
        case PARAM_ROOT_NOTE:
            info.id = PARAM_ROOT_NOTE;
            info.flags = CLAP_PARAM_IS_AUTOMATABLE | CLAP_PARAM_IS_STEPPED;
            info.min_value = ROOT_NOTE_MIN;
            info.max_value = ROOT_NOTE_MAX;
            info.default_value = ROOT_NOTE_DEFAULT;
            strncpy(info.name, "Root Note", CLAP_NAME_SIZE-1);
            strncpy(info.module, "", CLAP_PATH_SIZE-1);
            break;
        case PARAM_BANDWIDTH:
            info.id = PARAM_BANDWIDTH;
            info.flags = CLAP_PARAM_IS_AUTOMATABLE;
            info.min_value = BANDWIDTH_MIN;
            info.max_value = BANDWIDTH_MAX;
            info.default_value = BANDWIDTH_DEFAULT;
            strncpy(info.name, "Bandwidth", CLAP_NAME_SIZE-1);
            strncpy(info.module, "", CLAP_PATH_SIZE-1);
            break;
        case PARAM_ROOT_GAIN:
            info.id = PARAM_ROOT_GAIN;
            info.flags = CLAP_PARAM_IS_AUTOMATABLE;
            info.min_value = ROOT_GAIN_MIN;
            info.max_value = ROOT_GAIN_MAX;
            info.default_value = ROOT_GAIN_DEFAULT;
            strncpy(info.name, "Root Gain", CLAP_NAME_SIZE-1);
            strncpy(info.module, "", CLAP_PATH_SIZE-1);
            break;
        case PARAM_HARMONICS:
            info.id = PARAM_HARMONICS;
            info.flags = CLAP_PARAM_IS_AUTOMATABLE;
            info.min_value = HARMONICS_MIN;
            info.max_value = HARMONICS_MAX;
            info.default_value = HARMONICS_DEFAULT;
            strncpy(info.name, "Harmonics", CLAP_NAME_SIZE-1);
            strncpy(info.module, "", CLAP_PATH_SIZE-1);
            break;
        case PARAM_LO_OCTAVES:
            info.id = PARAM_LO_OCTAVES;
            info.flags = CLAP_PARAM_IS_AUTOMATABLE;
            info.min_value = LO_OCTAVES_MIN;
            info.max_value = LO_OCTAVES_MAX;
            info.default_value = LO_OCTAVES_DEFAULT;
            strncpy(info.name, "Low Octaves", CLAP_NAME_SIZE-1);
            strncpy(info.module, "", CLAP_PATH_SIZE-1);
            break;
        case PARAM_HI_OCTAVES:
            info.id = PARAM_HI_OCTAVES;
            info.flags = CLAP_PARAM_IS_AUTOMATABLE;
            info.min_value = HI_OCTAVES_MIN;
            info.max_value = HI_OCTAVES_MAX;
            info.default_value = HI_OCTAVES_DEFAULT;
            strncpy(info.name, "High Octaves", CLAP_NAME_SIZE-1);
            strncpy(info.module, "", CLAP_PATH_SIZE-1);
            break;
        case PARAM_SPREAD:
            info.id = PARAM_SPREAD;
            info.flags = CLAP_PARAM_IS_AUTOMATABLE;
            info.min_value = SPREAD_MIN;
            info.max_value = SPREAD_MAX;
            info.default_value = SPREAD_DEFAULT;
            strncpy(info.name, "Spread", CLAP_NAME_SIZE-1);
            strncpy(info.module, "", CLAP_PATH_SIZE-1);
            break;
        case PARAM_SATURATION:
            info.id = PARAM_SATURATION;
            info.flags = CLAP_PARAM_IS_AUTOMATABLE;
            info.min_value = SATURATION_MIN;
            info.max_value = SATURATION_MAX;
            info.default_value = SATURATION_DEFAULT;
            strncpy(info.name, "Saturation", CLAP_NAME_SIZE-1);
            strncpy(info.module, "", CLAP_PATH_SIZE-1);
            break;
        case PARAM_LPF_CUTOFF:
            info.id = PARAM_LPF_CUTOFF;
            info.flags = CLAP_PARAM_IS_AUTOMATABLE;
            info.min_value = LPF_CUTOFF_MIN;
            info.max_value = LPF_CUTOFF_MAX;
            info.default_value = LPF_CUTOFF_DEFAULT;
            strncpy(info.name, "LPF Cutoff", CLAP_NAME_SIZE-1);
            strncpy(info.module, "", CLAP_PATH_SIZE-1);
            break;
        case PARAM_WET_BOOST:
            info.id = PARAM_WET_BOOST;
            info.flags = CLAP_PARAM_IS_AUTOMATABLE;
            info.min_value = WET_BOOST_MIN;
            info.max_value = WET_BOOST_MAX;
            info.default_value = WET_BOOST_DEFAULT;
            strncpy(info.name, "Wet Boost", CLAP_NAME_SIZE-1);
            strncpy(info.module, "", CLAP_PATH_SIZE-1);
            break;
        case PARAM_DRY_WET:
            info.id = PARAM_DRY_WET;
            info.flags = CLAP_PARAM_IS_AUTOMATABLE;
            info.min_value = DRY_WET_MIN;
            info.max_value = DRY_WET_MAX;
            info.default_value = DRY_WET_DEFAULT;
            strncpy(info.name, "Dry/Wet", CLAP_NAME_SIZE-1);
            strncpy(info.module, "", CLAP_PATH_SIZE-1);
            break;
        default:
            return false;
    }
    // ensure null termination
    info.name[CLAP_NAME_SIZE-1] = '\0';
    info.module[CLAP_PATH_SIZE-1] = '\0';
    // copy to caller
    *param_info = info;
    return true;
}

static bool params_get_value(const clap_plugin *plugin, clap_id param_id, double *value) {
    minimal_plugin_t *p = (minimal_plugin_t*)plugin->plugin_data;
    
    switch (param_id) {
        case PARAM_FINE_TUNE:
            *value = p->fine_tune;
            return true;
        case PARAM_ROOT_NOTE:
            *value = p->root_note;
            return true;
        case PARAM_BANDWIDTH:
            *value = p->bandwidth;
            return true;
        case PARAM_ROOT_GAIN:
            *value = p->root_gain;
            return true;
        case PARAM_HARMONICS:
            *value = p->harmonics;
            return true;
        case PARAM_LO_OCTAVES:
            *value = p->lo_octaves;
            return true;
        case PARAM_HI_OCTAVES:
            *value = p->hi_octaves;
            return true;
        case PARAM_SPREAD:
            *value = p->spread;
            return true;
        case PARAM_SATURATION:
            *value = p->saturation;
            return true;
        case PARAM_LPF_CUTOFF:
            *value = p->lpf_cutoff;
            return true;
        case PARAM_WET_BOOST:
            *value = p->wet_boost;
            return true;
        case PARAM_DRY_WET:
            *value = p->dry_wet;
            return true;
    }
    return false;
}

static bool params_value_to_text(const clap_plugin *plugin, clap_id param_id, double value, char *display, uint32_t size) {
    switch (param_id) {
        case PARAM_FINE_TUNE:
            snprintf(display, size, "%+.1f cents", value);
            return true;
        case PARAM_ROOT_NOTE:
            {
                char note_name[16];
                get_note_name(value, note_name, sizeof(note_name));
                snprintf(display, size, "%s", note_name);
            }
            return true;
        case PARAM_BANDWIDTH:
            snprintf(display, size, "%.1f%%", value * 100.0);
            return true;
        case PARAM_ROOT_GAIN:
            snprintf(display, size, "%.0f%%", value * 100.0);
            return true;
        case PARAM_HARMONICS:
            snprintf(display, size, "%.0f%%", value * 100.0);
            return true;
        case PARAM_LO_OCTAVES:
            snprintf(display, size, "%.0f%%", value * 100.0);
            return true;
        case PARAM_HI_OCTAVES:
            snprintf(display, size, "%.0f%%", value * 100.0);
            return true;
        case PARAM_SPREAD:
            if (value == 0.0) {
                snprintf(display, size, "Center");
            } else if (value > 0.0) {
                snprintf(display, size, "R %.0f%%", value * 100.0);
            } else {
                snprintf(display, size, "L %.0f%%", (-value) * 100.0);
            }
            return true;
        case PARAM_SATURATION:
            if (value == 0.0) {
                snprintf(display, size, "Off");
            } else {
                snprintf(display, size, "%.1fx", value);
            }
            return true;
        case PARAM_LPF_CUTOFF:
            if (value >= 1000.0) {
                snprintf(display, size, "%.1f kHz", value / 1000.0);
            } else {
                snprintf(display, size, "%.0f Hz", value);
            }
            return true;
        case PARAM_WET_BOOST:
            if (value > 0.0) {
                snprintf(display, size, "+%.1f dB", value);
            } else {
                snprintf(display, size, "%.1f dB", value);
            }
            return true;
        case PARAM_DRY_WET:
            snprintf(display, size, "%.0f%%", value * 100.0);
            return true;
    }
    return false;
}

static bool params_text_to_value(const clap_plugin *plugin, clap_id param_id, const char *display, double *value) {
    return false;
}

static void params_flush(const clap_plugin *plugin, const clap_input_events_t *in, const clap_output_events_t *out) {
    minimal_plugin_t *p = (minimal_plugin_t*)plugin->plugin_data;
    
    uint32_t event_count = in->size(in);
    for (uint32_t i = 0; i < event_count; ++i) {
        const clap_event_header_t *header = in->get(in, i);
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
                    p->coefficients_need_update = true;
                    break;
                case PARAM_ROOT_GAIN:
                    p->root_gain = fmax(ROOT_GAIN_MIN, fmin(ROOT_GAIN_MAX, param_event->value));
                    p->coefficients_need_update = true;
                    break;
                case PARAM_HARMONICS:
                    p->harmonics = fmax(HARMONICS_MIN, fmin(HARMONICS_MAX, param_event->value));
                    p->coefficients_need_update = true;
                    break;
                case PARAM_LO_OCTAVES:
                    p->lo_octaves = fmax(LO_OCTAVES_MIN, fmin(LO_OCTAVES_MAX, param_event->value));
                    p->coefficients_need_update = true;
                    break;
                case PARAM_HI_OCTAVES:
                    p->hi_octaves = fmax(HI_OCTAVES_MIN, fmin(HI_OCTAVES_MAX, param_event->value));
                    p->coefficients_need_update = true;
                    break;
                case PARAM_SPREAD:
                    p->spread = fmax(SPREAD_MIN, fmin(SPREAD_MAX, param_event->value));
                    // Spread doesn't require coefficient updates
                    break;
                case PARAM_SATURATION:
                    p->saturation = fmax(SATURATION_MIN, fmin(SATURATION_MAX, param_event->value));
                    // Saturation doesn't require coefficient updates
                    break;
                case PARAM_LPF_CUTOFF:
                    p->lpf_cutoff = fmax(LPF_CUTOFF_MIN, fmin(LPF_CUTOFF_MAX, param_event->value));
                    p->coefficients_need_update = true;
                    break;
                case PARAM_WET_BOOST:
                    p->wet_boost = fmax(WET_BOOST_MIN, fmin(WET_BOOST_MAX, param_event->value));
                    // Wet boost doesn't require coefficient updates
                    break;
                case PARAM_DRY_WET:
                    p->dry_wet = fmax(DRY_WET_MIN, fmin(DRY_WET_MAX, param_event->value));
                    // Dry/wet doesn't require coefficient updates
                    break;
            }
        }
    }
}

const clap_plugin_params_t plugin_params = {
    params_count,
    params_get_info,
    params_get_value,
    params_value_to_text,
    params_text_to_value,
    params_flush
};

// audio ports extension
static uint32_t audio_ports_count(const clap_plugin *plugin, bool is_input) {
    return 1;
}

static bool audio_ports_get(const clap_plugin *plugin, uint32_t index, bool is_input, clap_audio_port_info_t *info) {
    if (index != 0) return false;
    
    info->id = is_input ? 0 : 1;
    strncpy(info->name, is_input ? "Audio Input" : "Audio Output", CLAP_NAME_SIZE - 1);
    info->name[CLAP_NAME_SIZE - 1] = '\0';
    info->flags = CLAP_AUDIO_PORT_IS_MAIN;
    info->channel_count = 2;
    info->port_type = CLAP_PORT_STEREO;
    info->in_place_pair = is_input ? CLAP_INVALID_ID : 0;
    
    return true;
}

const clap_plugin_audio_ports_t plugin_audio_ports = {
    audio_ports_count,
    audio_ports_get
};

// State extension

// Simple binary state format:
// - 4 bytes: magic number (0x45564853 = "EVHS")
// - 4 bytes: version number 
// - 12 x 8 bytes: parameter values as doubles

#define STATE_MAGIC 0x45564853  // "EVHS" in hex
#define STATE_VERSION 4  // Bumped for bug fix version

static bool state_save(const clap_plugin *plugin, const clap_ostream_t *stream) {
    if (!plugin || !stream) return false;
    
    minimal_plugin_t *p = (minimal_plugin_t*)plugin->plugin_data;
    
    // Write magic number and version
    uint32_t magic = STATE_MAGIC;
    uint32_t version = STATE_VERSION;
    
    if (stream->write(stream, &magic, sizeof(magic)) != sizeof(magic)) return false;
    if (stream->write(stream, &version, sizeof(version)) != sizeof(version)) return false;
    
    // Write all parameter values in order
    double params[PARAM_COUNT] = {
        p->fine_tune,
        p->root_note,
        p->bandwidth,
        p->root_gain,
        p->harmonics,
        p->lo_octaves,
        p->hi_octaves,
        p->spread,
        p->saturation,
        p->lpf_cutoff,
        p->wet_boost,
        p->dry_wet
    };
    
    size_t params_size = sizeof(params);
    if (stream->write(stream, params, params_size) != params_size) return false;
    
    return true;
}

static bool state_load(const clap_plugin *plugin, const clap_istream_t *stream) {
    if (!plugin || !stream) return false;
    
    minimal_plugin_t *p = (minimal_plugin_t*)plugin->plugin_data;
    
    // Read and verify magic number and version
    uint32_t magic, version;
    
    if (stream->read(stream, &magic, sizeof(magic)) != sizeof(magic)) return false;
    if (magic != STATE_MAGIC) return false;
    
    if (stream->read(stream, &version, sizeof(version)) != sizeof(version)) return false;
    if (version < 2 || version > STATE_VERSION) {
        // Handle version mismatch - for now just fail
        // You could add backward compatibility here
        return false;
    }
    
    // Read parameter values
    double params[PARAM_COUNT];
    size_t params_size = sizeof(params);
    if (stream->read(stream, params, params_size) != params_size) return false;
    
    // Assign values with bounds checking
    p->fine_tune = fmax(FINE_TUNE_MIN, fmin(FINE_TUNE_MAX, params[0]));
    p->root_note = fmax(ROOT_NOTE_MIN, fmin(ROOT_NOTE_MAX, params[1]));
    p->bandwidth = fmax(BANDWIDTH_MIN, fmin(BANDWIDTH_MAX, params[2]));
    p->root_gain = fmax(ROOT_GAIN_MIN, fmin(ROOT_GAIN_MAX, params[3]));
    p->harmonics = fmax(HARMONICS_MIN, fmin(HARMONICS_MAX, params[4]));
    p->lo_octaves = fmax(LO_OCTAVES_MIN, fmin(LO_OCTAVES_MAX, params[5]));
    p->hi_octaves = fmax(HI_OCTAVES_MIN, fmin(HI_OCTAVES_MAX, params[6]));
    p->spread = fmax(SPREAD_MIN, fmin(SPREAD_MAX, params[7]));
    p->saturation = fmax(SATURATION_MIN, fmin(SATURATION_MAX, params[8]));
    p->lpf_cutoff = fmax(LPF_CUTOFF_MIN, fmin(LPF_CUTOFF_MAX, params[9]));
    p->wet_boost = fmax(WET_BOOST_MIN, fmin(WET_BOOST_MAX, params[10]));
    p->dry_wet = fmax(DRY_WET_MIN, fmin(DRY_WET_MAX, params[11]));
    
    // Update smoothed parameter values to match loaded values
    p->bandwidth_smooth.current = p->bandwidth_smooth.target = p->bandwidth;
    p->root_gain_smooth.current = p->root_gain_smooth.target = p->root_gain;
    p->lo_octaves_smooth.current = p->lo_octaves_smooth.target = p->lo_octaves;
    p->hi_octaves_smooth.current = p->hi_octaves_smooth.target = p->hi_octaves;
    p->harmonics_smooth.current = p->harmonics_smooth.target = p->harmonics;
    p->saturation_smooth.current = p->saturation_smooth.target = p->saturation;
    p->lpf_cutoff_smooth.current = p->lpf_cutoff_smooth.target = p->lpf_cutoff;
    p->wet_boost_smooth.current = p->wet_boost_smooth.target = p->wet_boost;
    p->dry_wet_smooth.current = p->dry_wet_smooth.target = p->dry_wet;
    p->spread_smooth.current = p->spread_smooth.target = p->spread;
    
    // Reset smoothing flags
    p->bandwidth_smooth.active = false;
    p->root_gain_smooth.active = false;
    p->lo_octaves_smooth.active = false;
    p->hi_octaves_smooth.active = false;
    p->harmonics_smooth.active = false;
    p->saturation_smooth.active = false;
    p->lpf_cutoff_smooth.active = false;
    p->wet_boost_smooth.active = false;
    p->dry_wet_smooth.active = false;
    p->spread_smooth.active = false;
    p->any_smoothing_active = false;
    
    // Mark coefficients for update
    p->coefficients_need_update = true;
    
    return true;
}

const clap_plugin_state_t plugin_state = {
    state_save,
    state_load
};
