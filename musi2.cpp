/*
 * Copyright (c) 2015 Andrew Kelley
 *
 * This file is part of libsoundio, which is MIT licensed.
 * See http://opensource.org/licenses/MIT
 */
#include <soundio/soundio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <complex>
#include <iostream>

using namespace std;


static int usage(char *exe) {
    fprintf(stderr, "Usage: %s [options]\n"
            "Options:\n"
            "  [--backend dummy|alsa|pulseaudio|jack|coreaudio|wasapi]\n"
            "  [--device id]\n"
            "  [--raw]\n"
            "  [--name stream_name]\n"
            "  [--latency seconds]\n"
            "  [--sample-rate hz]\n"
            "  [--nMelodies >0]\n"
            "  [--nDrums >0]\n"
            "  [--probDoDrum 0-100]\n"
            "  [--probDoKey 0-100]\n"
            "  [--probDoNextKey 0-100]\n"
            , exe);
    return 1;
}


int bufferIndex = 0;
signed short int buffer[18 * 44100];

int probDoDrum = 50;
int probDoKey = 80;
int probDoNextKey = 80;


int tonesXY[20];
int numTonesXY = 0;

int scale[12];


class Instrument {
	public:
	float leftRight;
	float sr;
	float a;
	bool b;
	float c;
	bool d;
	
	Instrument(int srate,float f,float g,bool x) {
		leftRight = (float)rand() / (float)RAND_MAX;
		sr = srate;
		a = f * rand() / RAND_MAX;
		b = rand() % 2 == 0;
		c = g * rand() / RAND_MAX;
		if(x)
			d = rand() % 2 == 0;
		else
			d = false;
	}
	
	virtual int render(int index,int sz,int scaleIndex,double fr,float vol);
};
int Instrument::render(int index,int sz,int scaleIndex,double fr,float vol) {
		/*a = 0.29 * rand() / RAND_MAX;
		b = true;
		c = 0.09 * rand() / RAND_MAX;
		d = true;*/

	int n = rand() % 500;
	
	if(scaleIndex >= 0) {
	////////////////////////////////////////////////
	// Harmonisierung 
	// Asure that there is no 1-half-tone-interval at all
	const double xy = pow(2.0,1.0 / 12.0);
	const double xy2 = pow(2.0,0.8 / 12.0);
	const double xy3 = pow(2.0,1.2 / 12.0);
	
	bool found;
	do {
		found = false;
		fr = pow(2.0,scale[scaleIndex] / 12.0) * sr / 44100.0;	
		
		for(int i = 0;i < numTonesXY;++i) {
		fr = pow(2.0,scale[scaleIndex] / 12.0) * sr / 44100.0;	
			double f = fr / scale[tonesXY[i]];
			//printf("%f %f %f\n",fr,xy2,xy3);
			if((f > xy2 && f < xy3) || (1 / f > xy2 && 1 / f < xy3)) {
				
				if(rand() % 2)
					scaleIndex += 1;
				else
					scaleIndex -= 1;
				
				if(scaleIndex < 0) scaleIndex = 0;
				if(scaleIndex >= sizeof(scale) / sizeof(scale[0])) scaleIndex = sizeof(scale) / sizeof(scale[0]) - 1;
				
				i = 0;
				found = true;
			} 
		}
	}
	while(found);
	
	tonesXY[numTonesXY] = scaleIndex;
	numTonesXY += 1;
	//
	/////////////////////////////////////////////////
	} 
	
	for(int i = 0;i < sz;++i) 
		if(i + index + n < sizeof(buffer) / sizeof(buffer[0]))
			{
				double samp = 0.1 * vol * 12000.0 * sin(i * fr * 440.0 * M_PI / 44100.0);
				if(b) samp *= cos(a * i * fr * 440.0 * M_PI / 44100.0) / 2 + 0.5;
				if(c) samp *= cos(d * i * 440.0 * M_PI / 44100.0) / 2 + 0.5;
				buffer[i + index + n] += samp * leftRight;
				buffer[i + index + n] += samp * (1 - leftRight);
			}
	
	return (int)(20000.0 * fr);
}




struct Drum {
	Instrument ins;
	
	Drum() : ins(44100,0.2,0.1,true) {
	}
	
	void render(int index) {
		if(rand() % 100 < probDoDrum)
			ins.render(index,6000,-1,0.75,8.0);
	}
};




void makeScale() {
	int sc[] = {0,2,4,5,7,9,11};
	/*int l = sizeof(sc) / sizeof(sc[0]);
	int sc0[l];
	
	for(int i = 0;i < l - 1;++i) sc0[i] = sc[i + 1] - sc[i];
	sc0[l - 1] = 12 + sc[0] - sc[l - 1];
	int m = rand() % 12;
	
	for(int i = 0;i < m;++i) {
		int swp = sc0[0];
		for(int j = 1;j < l;++j) sc0[j - 1] = sc0[j];
		sc0[l - 1] = swp;
	}
	for(int j = 1;j < l;++j) sc[j] = sc0[j - 1] + sc[j - 1];
	printf("[");
	for(int i = 0;i < sizeof(sc0) / sizeof(sc0[0]);++i) printf("%d,",sc0[i]);
	printf("]\n");*/
	
	
	int n = rand() % 8;
	
	int k = 0;
	for(int j = 0;j < 5;++j)
	for(int i = 0;i < sizeof(sc) / sizeof(sc[0]);++i) 
	if(k < sizeof(scale) / sizeof(scale[0])) {
		scale[k] = sc[i] + 12 * j + n;
		k += 1;
	}
	
	printf("[");
	for(int i = 0;i < sizeof(scale) / sizeof(scale[0]);++i) printf("%d,",scale[i]);
	printf("]\n");
}


struct Melody {
	Instrument *ins;
	float tones[32];
	int index;
	float min,max;
	
	
	void fractal(int start,int end,double decay,double amp) {
		int m = (start + end) / 2;
		if(tones[m] == 0) {
			tones[m] = (tones[start] + tones[end]) / 2.0 + ((float)rand() / (float)RAND_MAX) * amp;
			fractal(start,m,decay,-decay * amp);
			fractal(m,end,decay,-decay * amp);
		}
	}
	
	
	
	void init() {
		index = 0;
		for(int i = 0;i < sizeof(tones) / sizeof(tones[0]);++i) tones[i] = 0;	
		tones[0] = (sizeof(scale) / sizeof(scale[0])) * (float)rand() / (float)RAND_MAX;
		tones[sizeof(tones) / sizeof(tones[0]) - 1] = (sizeof(scale) / sizeof(scale[0])) * (float)rand() / (float)RAND_MAX;
		for(int i = 0;i < sizeof(tones) / sizeof(tones[0]);++i) printf("%f,",tones[i]);	
		printf("\n");	
			
		fractal(0,sizeof(tones) / sizeof(tones[0]) - 1,0.99,sizeof(scale) / sizeof(scale[0]) * (rand() % 2 ? 1.0 : -1.0));
		for(int i = 0;i < sizeof(tones) / sizeof(tones[0]);++i) printf("%f,",tones[i]);	
		printf("\n");	
		
		min = 100000;
		max = -100000;
		for(int i = 0;i < sizeof(tones) / sizeof(tones[0]);++i) if(tones[i] > max) max = tones[i];		
		for(int i = 0;i < sizeof(tones) / sizeof(tones[0]);++i) if(tones[i] < min) min = tones[i];	
		for(int i = 0;i < sizeof(tones) / sizeof(tones[0]);++i) printf("%f,",tones[i]);	
		printf("\n");	
		
		printf("%f\n",min);	
		printf("%f\n",max);	
	}
	
	Melody() : ins(new Instrument(44100,0.0,0.00000,true)) {
		init();
	}
	
	Melody(Instrument *in) : ins(in) {
		init();
	}
	
	void render(int ind) {
		if(rand() % 100 < probDoKey) {
			int k = (int)((tones[index] - min) / (max - min) * (sizeof(scale) / sizeof(scale[0]) - 1));
			printf("%d  ",k);
			
			ins->render(ind,20000,k,0.0,1.5);
			
			if(rand() % 100 < probDoNextKey) {
				index += 1;
				
				if(index >= sizeof(tones) / sizeof(tones[0]))
					index = 0;
			}
		}
	}
};



static void write_sample_s16ne(char *ptr, double sample) {
    int16_t *buf = (int16_t *)ptr;
    double range = (double)INT16_MAX - (double)INT16_MIN;
    double val = sample * range / 2.0;
    *buf = val;
}
static void write_sample_s32ne(char *ptr, double sample) {
    int32_t *buf = (int32_t *)ptr;
    double range = (double)INT32_MAX - (double)INT32_MIN;
    double val = sample * range / 2.0;
    *buf = val;
}
static void write_sample_float32ne(char *ptr, double sample) {
    float *buf = (float *)ptr;
    *buf = sample;
}
static void write_sample_float64ne(char *ptr, double sample) {
    double *buf = (double *)ptr;
    *buf = sample;
}





bool doRender = true;
int nDrums = 4;
int nMelodies = 4;



static void (*write_sample)(char *ptr, double sample);
static const double PI = 3.14159265358979323846264338328;
static double seconds_offset = 0.0;
static volatile bool want_pause = false;
static void write_callback(struct SoundIoOutStream *outstream, int frame_count_min, int frame_count_max) {
    double float_sample_rate = outstream->sample_rate;
    double seconds_per_frame = 1.0 / float_sample_rate;
    struct SoundIoChannelArea *areas;
    int err;
    int frames_left = frame_count_max;
    for (;;) {
        int frame_count = frames_left;
        if ((err = soundio_outstream_begin_write(outstream, &areas, &frame_count))) {
            fprintf(stderr, "unrecoverable stream error: %s\n", soundio_strerror(err));
            exit(1);
        }
        if (!frame_count)
            break;
        const struct SoundIoChannelLayout *layout = &outstream->layout;
        double pitch = 440.0;
        double radians_per_second = pitch * 2.0 * PI;
        for (int frame = 0; frame < frame_count; frame += 1) {
            double sample = sin((seconds_offset + frame * seconds_per_frame) * radians_per_second);
            sample = (double)buffer[bufferIndex] / 32768.0;
            
            bufferIndex += 1;
            if(doRender || bufferIndex >= sizeof(buffer) / sizeof(buffer[0])) {
            	doRender = false;
            	
            	bufferIndex = 0; 
            	
	    	Drum ds[nDrums];
	    	makeScale();
	    	Melody ms[nMelodies];
	    	
	    	printf("rendering\n");
            for(int k = 0;k < sizeof(buffer) / sizeof(buffer[0]);++k) buffer[k] = 0;

			int cnt = 0;
    		for(int i = 0;i < sizeof(buffer) / sizeof(buffer[0]);i += 5500) {
    			numTonesXY = 0;

    			for(int j = 0;j < sizeof(ds) / sizeof(ds[0]);++j) ds[j].render(i);
    			
    			numTonesXY = 0;
    			
    			if(cnt == 0)
    			for(int j = 0;j < sizeof(ms) / sizeof(ms[0]);++j) ms[j].render(i);
    			
    			
    			cnt += 1;
    			if(cnt >= 2) cnt = 0;
    		}
	    	printf("playing\n");
            }
            
            for (int channel = 0; channel < layout->channel_count; channel += 1) {
                write_sample(areas[channel].ptr, sample);
                areas[channel].ptr += areas[channel].step;
            }
        }
        seconds_offset = fmod(seconds_offset + seconds_per_frame * frame_count, 1.0);
        if ((err = soundio_outstream_end_write(outstream))) {
            if (err == SoundIoErrorUnderflow)
                return;
            fprintf(stderr, "unrecoverable stream error: %s\n", soundio_strerror(err));
            exit(1);
        }
        frames_left -= frame_count;
        if (frames_left <= 0)
            break;
    }
    soundio_outstream_pause(outstream, want_pause);
}
static void underflow_callback(struct SoundIoOutStream *outstream) {
    static int count = 0;
    fprintf(stderr, "underflow %d\n", count++);
}

void FFT(complex<double> *data, unsigned N){
	unsigned butterflySize;   // size for actual butterfly calculation
	int i, j, k;              // local index variables
	complex<double> wActual;  // actual rotation factor
	complex<double> wStep;    // step rotation factors
	complex<double> tmp;      // temp. value for butterfly calculation
	// loop over all level of FFT
	for(butterflySize=N/2; butterflySize>0; butterflySize/=2) {
		// evaluate angle step and set first angle
		wStep = complex<double>(cos(-M_PI/butterflySize),sin(-M_PI/butterflySize));
		wActual = complex<double>(1, 0);
		// loop over number of butterflys
		for(j=0; j<butterflySize; j++) {
			// loop over number of FFTs
			for(i=j; i<N; i+=2*butterflySize) {
				// get index of second element
				k = i+butterflySize;
				// perform butterfly calculation
				tmp = data[i];          // store one element
				data[i] += data[k];     // take sum
				data[k] = tmp-data[k];  // take difference
				data[k] *= wActual;     // multiply with rotation factor
			}
			// evaluate next rotation factor
			wActual *= wStep;
		}
	}
	
	// perform bit reversal
	j = 0;
	for(i=0; i<N; i++) {
		if(j>i) {
			// swap numbers
			tmp = data[i];
			data[i] = data[j];
			data[j] = tmp;
		}
		k = N/2;
		while(k>=2 && j>=k) {
			j -= k;
			k /= 2;
		}
		j += k;
	}
}

void RFFT(complex<double> *data, unsigned N){
	unsigned butterflySize;   // size for actual butterfly calculation
	int i, j, k;              // local index variables
	complex<double> wActual;  // actual rotation factor
	complex<double> wStep;    // step rotation factors
	complex<double> tmp;      // temp. value for butterfly calculation
	// loop over all level of FFT
	for(butterflySize=N/2; butterflySize>0; butterflySize/=2) {
		// evaluate angle step and set first angle
		wStep = complex<double>(cos(+M_PI/butterflySize),sin(+M_PI/butterflySize));
		wActual = complex<double>(1, 0);
		// loop over number of butterflys
		for(j=0; j<butterflySize; j++) {
			// loop over number of FFTs
			for(i=j; i<N; i+=2*butterflySize) {
				// get index of second element
				k = i+butterflySize;
				// perform butterfly calculation
				tmp = data[i];          // store one element
				data[i] += data[k];     // take sum
				data[k] = tmp-data[k];  // take difference
				data[k] *= wActual;     // multiply with rotation factor
			}
			// evaluate next rotation factor
			wActual *= wStep;
		}
	}
	
	// perform bit reversal
	j = 0;
	for(i=0; i<N; i++) {
		if(j>i) {
			// swap numbers
			tmp = data[i];
			data[i] = data[j];
			data[j] = tmp;
		}
		k = N/2;
		while(k>=2 && j>=k) {
			j -= k;
			k /= 2;
		}
		j += k;
	}
}

int main(int argc, char **argv) {
	srand(time(0));

	
    char *exe = argv[0];
    enum SoundIoBackend backend = SoundIoBackendNone;
    char *device_id = NULL;
    bool raw = false;
    char *stream_name = NULL;
    double latency = 0.0;
    int sample_rate = 44100;
    for (int i = 1; i < argc; i += 1) {
        char *arg = argv[i];
        if (arg[0] == '-' && arg[1] == '-') {
            if (strcmp(arg, "--raw") == 0) {
                raw = true;
            } else {
                i += 1;
                if (i >= argc) {
                    return usage(exe);
                } else if (strcmp(arg, "--backend") == 0) {
                    if (strcmp(argv[i], "dummy") == 0) {
                        backend = SoundIoBackendDummy;
                    } else if (strcmp(argv[i], "alsa") == 0) {
                        backend = SoundIoBackendAlsa;
                    } else if (strcmp(argv[i], "pulseaudio") == 0) {
                        backend = SoundIoBackendPulseAudio;
                    } else if (strcmp(argv[i], "jack") == 0) {
                        backend = SoundIoBackendJack;
                    } else if (strcmp(argv[i], "coreaudio") == 0) {
                        backend = SoundIoBackendCoreAudio;
                    } else if (strcmp(argv[i], "wasapi") == 0) {
                        backend = SoundIoBackendWasapi;
                    } else {
                        fprintf(stderr, "Invalid backend: %s\n", argv[i]);
                        return 1;
                    }
                } else if (strcmp(arg, "--device") == 0) {
                    device_id = argv[i];
                } else if (strcmp(arg, "--probDoDrum") == 0) {
                    probDoDrum = atoi(argv[i]);
                } else if (strcmp(arg, "--probDoKey") == 0) {
                    probDoKey = atoi(argv[i]);
                } else if (strcmp(arg, "--probDoNextKey") == 0) {
                    probDoNextKey = atoi(argv[i]);
                } else if (strcmp(arg, "--nDrums") == 0) {
                    nDrums = atoi(argv[i]);
                } else if (strcmp(arg, "--nMelodies") == 0) {
                    nMelodies = atoi(argv[i]);
                } else if (strcmp(arg, "--name") == 0) {
                    stream_name = argv[i];
                } else if (strcmp(arg, "--latency") == 0) {
                    latency = atof(argv[i]);
                } else if (strcmp(arg, "--sample-rate") == 0) {
                    sample_rate = atoi(argv[i]);
                } else {
                    return usage(exe);
                }
            }
        } else {
            return usage(exe);
        }
    }
    struct SoundIo *soundio = soundio_create();
    if (!soundio) {
        fprintf(stderr, "out of memory\n");
        return 1;
    }
    int err = (backend == SoundIoBackendNone) ?
        soundio_connect(soundio) : soundio_connect_backend(soundio, backend);
    if (err) {
        fprintf(stderr, "Unable to connect to backend: %s\n", soundio_strerror(err));
        return 1;
    }
    fprintf(stderr, "Backend: %s\n", soundio_backend_name(soundio->current_backend));
    soundio_flush_events(soundio);
    int selected_device_index = -1;
    if (device_id) {
        int device_count = soundio_output_device_count(soundio);
        for (int i = 0; i < device_count; i += 1) {
            struct SoundIoDevice *device = soundio_get_output_device(soundio, i);
            bool select_this_one = strcmp(device->id, device_id) == 0 && device->is_raw == raw;
            soundio_device_unref(device);
            if (select_this_one) {
                selected_device_index = i;
                break;
            }
        }
    } else {
        selected_device_index = soundio_default_output_device_index(soundio);
    }
    if (selected_device_index < 0) {
        fprintf(stderr, "Output device not found\n");
        return 1;
    }
    struct SoundIoDevice *device = soundio_get_output_device(soundio, selected_device_index);
    if (!device) {
        fprintf(stderr, "out of memory\n");
        return 1;
    }
    fprintf(stderr, "Output device: %s\n", device->name);
    if (device->probe_error) {
        fprintf(stderr, "Cannot probe device: %s\n", soundio_strerror(device->probe_error));
        return 1;
    }
    struct SoundIoOutStream *outstream = soundio_outstream_create(device);
    if (!outstream) {
        fprintf(stderr, "out of memory\n");
        return 1;
    }
    outstream->write_callback = write_callback;
    outstream->underflow_callback = underflow_callback;
    outstream->name = stream_name;
    outstream->software_latency = latency;
    outstream->sample_rate = sample_rate;
    if (soundio_device_supports_format(device, SoundIoFormatFloat32NE)) {
        outstream->format = SoundIoFormatFloat32NE;
        write_sample = write_sample_float32ne;
    } else if (soundio_device_supports_format(device, SoundIoFormatFloat64NE)) {
        outstream->format = SoundIoFormatFloat64NE;
        write_sample = write_sample_float64ne;
    } else if (soundio_device_supports_format(device, SoundIoFormatS32NE)) {
        outstream->format = SoundIoFormatS32NE;
        write_sample = write_sample_s32ne;
    } else if (soundio_device_supports_format(device, SoundIoFormatS16NE)) {
        outstream->format = SoundIoFormatS16NE;
        write_sample = write_sample_s16ne;
    } else {
        fprintf(stderr, "No suitable device format available.\n");
        return 1;
    }
    if ((err = soundio_outstream_open(outstream))) {
        fprintf(stderr, "unable to open device: %s", soundio_strerror(err));
        return 1;
    }
    fprintf(stderr, "Software latency: %f\n", outstream->software_latency);
    fprintf(stderr,
            "'p\\n' - pause\n"
            "'u\\n' - unpause\n"
            "'P\\n' - pause from within callback\n"
            "'c\\n' - clear buffer\n"
            "'q\\n' - quit\n");
    if (outstream->layout_error)
        fprintf(stderr, "unable to set channel layout: %s\n", soundio_strerror(outstream->layout_error));
    if ((err = soundio_outstream_start(outstream))) {
        fprintf(stderr, "unable to start device: %s\n", soundio_strerror(err));
        return 1;
    }
    
    bufferIndex = 0;
    
    for (;;) {
    	
        soundio_flush_events(soundio);
        int c = getc(stdin);
        if (c == 'p') {
            fprintf(stderr, "pausing result: %s\n",
                    soundio_strerror(soundio_outstream_pause(outstream, true)));
        } else if (c == 'P') {
            want_pause = true;
        } else if (c == 'u') {
            want_pause = false;
            fprintf(stderr, "unpausing result: %s\n",
                    soundio_strerror(soundio_outstream_pause(outstream, false)));
        } else if (c == 'c') {
            fprintf(stderr, "clear buffer result: %s\n",
                    soundio_strerror(soundio_outstream_clear_buffer(outstream)));
        } else if (c == 'q') {
            break;
        } else if (c == '\r' || c == '\n') {
            // ignore
        } else {
            fprintf(stderr, "Unrecognized command: %c\n", c);
        }
    }
    soundio_outstream_destroy(outstream);
    soundio_device_unref(device);
    soundio_destroy(soundio);
    return 0;
}
