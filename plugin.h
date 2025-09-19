#pragma once
#include <clap/clap.h>

// plugin identification
#define PLUGIN_ID      "com.evan.harmonicisolator"
#define PLUGIN_NAME    "Evan's Harmonic Isolator"  
#define PLUGIN_VENDOR  "Evan"
#define PLUGIN_VERSION "0.9.7"  // Elliptic filter release
#define PLUGIN_DESC    "A tunable harmonic isolator with stepped frequency selection - elliptic AA filters"

// parameter definitions - FINE_TUNE is now at position 0!
enum {
    PARAM_FINE_TUNE = 0,  // NEW! Placeholder parameter
    PARAM_ROOT_NOTE,      // now 1
    PARAM_BANDWIDTH,      // now 2
    PARAM_ROOT_GAIN,      // now 3
    PARAM_HARMONICS,      // now 4
    PARAM_LO_OCTAVES,     // now 5
    PARAM_HI_OCTAVES,     // now 6
    PARAM_SPREAD,         // now 7
    PARAM_SATURATION,     // now 8
    PARAM_LPF_CUTOFF,     // now 9
    PARAM_WET_BOOST,      // now 10
    PARAM_DRY_WET,        // now 11
    PARAM_COUNT           // should be 12
};

// Fine tune parameter (cents)
#define FINE_TUNE_MIN -100.0
#define FINE_TUNE_MAX 100.0
#define FINE_TUNE_DEFAULT 0.0

#define ROOT_NOTE_MIN 0.0
#define ROOT_NOTE_MAX 39.0  // 40 semitones from A3 to C6
#define ROOT_NOTE_DEFAULT 12.0  // A4 (440 Hz)

#define BANDWIDTH_MIN 0.01
#define BANDWIDTH_MAX 1.0
#define BANDWIDTH_DEFAULT 0.4  // 40%

#define ROOT_GAIN_MIN 0.0
#define ROOT_GAIN_MAX 1.0
#define ROOT_GAIN_DEFAULT 0.6  // 60%

#define LO_OCTAVES_MIN 0.0
#define LO_OCTAVES_MAX 1.0
#define LO_OCTAVES_DEFAULT 0.1

#define HI_OCTAVES_MIN 0.0
#define HI_OCTAVES_MAX 1.0
#define HI_OCTAVES_DEFAULT 0.2  // 20% - changed from 30%

#define HARMONICS_MIN 0.0
#define HARMONICS_MAX 1.0
#define HARMONICS_DEFAULT 0.1  // 10%

#define SATURATION_MIN 0.0
#define SATURATION_MAX 10.0
#define SATURATION_DEFAULT 2.4  // subtle saturation by default

#define LPF_CUTOFF_MIN 40.0
#define LPF_CUTOFF_MAX 20000.0
#define LPF_CUTOFF_DEFAULT 8000.0  // default fairly open

#define DRY_WET_MIN 0.0
#define DRY_WET_MAX 1.0
#define DRY_WET_DEFAULT 1.0  // default to 100% wet

#define WET_BOOST_MIN -10.0
#define WET_BOOST_MAX 20.0
#define WET_BOOST_DEFAULT 6.0

#define SPREAD_MIN -1.0
#define SPREAD_MAX 1.0
#define SPREAD_DEFAULT 0.0  // center by default

// Oversampling constants
#define OVERSAMPLE_FACTOR 2
#define AA_FILTER_ORDER 8  // 8th order elliptic filter for anti-aliasing

// Math constants
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Optimization constants
#define MAX_ACTIVE_FILTERS 16  // Safety limit for active filters

// Parameter smoothing constants
#define PARAM_SMOOTH_TIME_MS 10.0f  // 10ms smoothing time

// Multiple bandpass filters - one per semitone in octave
#define NUM_SEMITONES 12
#define NUM_OCTAVES 12  // Increased from 8 to handle all possible filters

// Parameter smoothing structure
typedef struct {
    float current;     // current smoothed value
    float target;      // target value to reach
    bool active;       // is this parameter currently smoothing?
} param_smooth_t;

// Anti-aliasing filter structure (8th order = 4 biquads)
typedef struct {
    float x1[4], x2[4];  // input delay lines for 4 biquads
    float y1[4], y2[4];  // output delay lines for 4 biquads
    float b0[4], b1[4], b2[4];  // numerator coefficients
    float a1[4], a2[4];         // denominator coefficients (a0 = 1)
} aa_filter_t;

// Active filter tracking
typedef struct {
    int octave_idx;
    int semitone_idx;
} active_filter_t;

// plugin state structure
typedef struct {
    clap_host_t *host;
    
    // parameter values (current)
    double fine_tune;     // NEW! fine tune in cents
    double root_note;     // semitone offset from A3
    double bandwidth;
    double root_gain;     // renamed from root_level
    double lo_octaves;    // split from octaves
    double hi_octaves;    // split from octaves
    double harmonics;
    double saturation;    // NEW! saturation drive amount
    double lpf_cutoff;
    double wet_boost;
    double dry_wet;
    double spread;        // NEW! stereo spread amount (bipolar)
    
    // Parameter smoothing system
    param_smooth_t bandwidth_smooth;
    param_smooth_t root_gain_smooth;
    param_smooth_t lo_octaves_smooth;
    param_smooth_t hi_octaves_smooth;
    param_smooth_t harmonics_smooth;
    param_smooth_t saturation_smooth;
    param_smooth_t lpf_cutoff_smooth;
    param_smooth_t wet_boost_smooth;
    param_smooth_t dry_wet_smooth;
    param_smooth_t spread_smooth;
    
    float smooth_coeff;              // calculated from sample rate
    bool any_smoothing_active;       // quick check flag
    
    // DSP state
    double sample_rate;
    double oversampled_rate;  // sample_rate * OVERSAMPLE_FACTOR
    uint32_t max_block_size;  // from activate() for buffer allocation
    
    // Active filter tracking - optimized for better cache performance
    active_filter_t active_filters[MAX_ACTIVE_FILTERS];
    int num_active_filters;
    
    // Anti-aliasing filters
    aa_filter_t upsample_aa_L;    // anti-aliasing for upsampling (L)
    aa_filter_t upsample_aa_R;    // anti-aliasing for upsampling (R)
    aa_filter_t downsample_aa_L;  // anti-aliasing for downsampling (L)
    aa_filter_t downsample_aa_R;  // anti-aliasing for downsampling (R)
    
    // Oversampled processing buffers - SIMD aligned
    float *oversample_buffer_L;
    float *oversample_buffer_R;
    size_t oversample_buffer_size;
    
    // LPF state (adaptive multi-pole filters for each channel) - now at oversampled rate
    float lpf_state_L;
    float lpf_state_R;
    float lpf2_state_L;  // second pole state
    float lpf2_state_R;  // second pole state
    float lpf3_state_L;  // third pole state
    float lpf3_state_R;  // third pole state
    float lpf_coeff;
    
    // Biquad coefficients for each filter - optimized layout
    float bp_b0[NUM_OCTAVES][NUM_SEMITONES];
    float bp_b1[NUM_OCTAVES][NUM_SEMITONES];
    float bp_b2[NUM_OCTAVES][NUM_SEMITONES];
    float bp_a1[NUM_OCTAVES][NUM_SEMITONES];
    float bp_a2[NUM_OCTAVES][NUM_SEMITONES];
    
    // Filter state for left channel
    float bp_x1_L[NUM_OCTAVES][NUM_SEMITONES];
    float bp_x2_L[NUM_OCTAVES][NUM_SEMITONES];
    float bp_y1_L[NUM_OCTAVES][NUM_SEMITONES];
    float bp_y2_L[NUM_OCTAVES][NUM_SEMITONES];
    
    // Filter state for right channel  
    float bp_x1_R[NUM_OCTAVES][NUM_SEMITONES];
    float bp_x2_R[NUM_OCTAVES][NUM_SEMITONES];
    float bp_y1_R[NUM_OCTAVES][NUM_SEMITONES];
    float bp_y2_R[NUM_OCTAVES][NUM_SEMITONES];
    
    // parameter change tracking
    double last_fine_tune;
    double last_root_note;
    double last_bandwidth;
    double last_lo_octaves;
    double last_hi_octaves;
    double last_harmonics;
    double last_lpf_cutoff;
    double last_root_gain;  // ADDED
    bool coefficients_need_update;
    
} minimal_plugin_t;

#ifdef __cplusplus
extern "C" {
#endif

// plugin lifecycle functions (plugin.cpp)
bool minimal_init(const clap_plugin *plugin);
void minimal_destroy(const clap_plugin *plugin);
bool minimal_activate(const clap_plugin *plugin, double sample_rate, uint32_t min_frames, uint32_t max_frames);
void minimal_deactivate(const clap_plugin *plugin);
bool minimal_start_processing(const clap_plugin *plugin);
void minimal_stop_processing(const clap_plugin *plugin);
void minimal_reset(const clap_plugin *plugin);
clap_process_status minimal_process(const clap_plugin *plugin, const clap_process_t *process);
const void *minimal_get_extension(const clap_plugin *plugin, const char *id);
void minimal_on_main_thread(const clap_plugin *plugin);

// DSP functions (dsp.cpp)
void update_filter_coefficients(minimal_plugin_t *p);
void reset_all_filter_states(minimal_plugin_t *p);
void init_oversampling(minimal_plugin_t *p);
void cleanup_oversampling(minimal_plugin_t *p);
void reset_aa_filters(minimal_plugin_t *p);
float process_aa_filter(float input, aa_filter_t *filter);
void update_parameter_smoothing(minimal_plugin_t *p);
void trigger_parameter_smoothing(minimal_plugin_t *p);
void process_parameter_smoothing(minimal_plugin_t *p, uint32_t frames);

// Filter processing functions (filters.cpp)
void calc_constant_power_pan(float spread, float *pan_L, float *pan_R);
float process_adaptive_lpf(float combined_output, float &lpf_state, float &lpf2_state, float &lpf3_state, float lpf_coeff, float lpf_cutoff_smooth);
void upsample_2x(const float *input, float *output, uint32_t input_frames, aa_filter_t *aa_filter);
void downsample_2x(const float *input, float *output, uint32_t input_frames, aa_filter_t *aa_filter);
void process_channel_filters(minimal_plugin_t *p, float *buffer, uint32_t frames, bool is_left_channel);

// Parameter handling functions (params.cpp)
extern const clap_plugin_params_t plugin_params;
extern const clap_plugin_audio_ports_t plugin_audio_ports;
extern const clap_plugin_state_t plugin_state;

// Utility functions (utils.cpp)
double semitone_to_freq(double semitone);
void get_note_name(double semitone, char* buffer, size_t size);

#ifdef __cplusplus
}
#endif
