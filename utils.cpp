#include "plugin.h"
#include <cmath>
#include <cstdio>

// Note names for display
static const char* note_names[] = {
    "A", "A#", "B", "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#"
};

// Helper function to convert semitone offset from A3 to frequency
double semitone_to_freq(double semitone) {
    // A3 = 220 Hz, semitone is offset from A3
    return 220.0 * pow(2.0, semitone / 12.0);
}

// Helper function to get note name and octave from semitone offset
void get_note_name(double semitone, char* buffer, size_t size) {
    int semi = (int)(semitone + 0.5);  // round to nearest
    int note_in_octave = semi % 12;
    int octave = 3 + (semi + 3) / 12;  // +3 because we start at A (9th note)
    
    snprintf(buffer, size, "%s%d", note_names[note_in_octave], octave);
}
