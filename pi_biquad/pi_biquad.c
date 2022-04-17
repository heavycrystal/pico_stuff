#include        <stdint.h>
#include        <stdio.h>
#include        <string.h>

#include        "hardware/vreg.h"
#include        "pico/stdio_usb.h"
#include        "pico/stdlib.h"
#include	    "pico/multicore.h"
#include        "pico/time.h"
#include        "pico/types.h"

#define         LED_PIN         PICO_DEFAULT_LED_PIN
#define         ARRAY_SIZE      32768U
#define         ITERATIONS      4096U

const int16_t a0 = 16384;
const int16_t a1 = -32768;
const int16_t a2 = 16384;
const int16_t b1 = -25576;
const int16_t b2 = 10508;

void core1_main(void)
{
    absolute_time_t start_time;
    size_t duration_ms;
    volatile int16_t in[ARRAY_SIZE];
    volatile int16_t out[ARRAY_SIZE];
    int16_t z1, z2;
    int16_t outTemp;
    int16_t inTemp;
    uint64_t counter = 1;

    while(1) 
    {
        start_time = get_absolute_time();
        for(size_t loop_var = 0u; loop_var < ITERATIONS; loop_var++)
        {
            for(size_t i = 0; i < ARRAY_SIZE; i++)
            {
                inTemp = in[i];
                outTemp = inTemp * a0 + z1;
                z1 = inTemp * a1 + z2 - b1 * outTemp;
                z2 = inTemp * a2 - b2 * outTemp;
                out[i] = outTemp;
            }
        }
        duration_ms = absolute_time_diff_us(start_time, get_absolute_time()) / 1000;
        printf("[Core #1] Finished iteration #%llu.\n", counter);
        printf("[Core #1] Finished %u samples in %zu milliseconds.\n", ARRAY_SIZE * ITERATIONS, duration_ms);
        printf("[Core #1] %zu kiloSamples per second.\n", ARRAY_SIZE * ITERATIONS / duration_ms);
        counter = counter + 1;
    }
}

int main(void)
{
    gpio_init(PICO_DEFAULT_LED_PIN);
    gpio_set_dir(LED_PIN, GPIO_OUT);
    gpio_put(LED_PIN, 1);

    stdio_usb_init();
    multicore_launch_core1(core1_main);

    absolute_time_t start_time;
    size_t duration_ms;
    volatile int16_t in[ARRAY_SIZE];
    volatile int16_t out[ARRAY_SIZE];
    int16_t z1, z2;
    int16_t outTemp;
    int16_t inTemp;
    uint64_t counter = 1;

    while(1) 
    {
        start_time = get_absolute_time();
        for(size_t loop_var = 0u; loop_var < ITERATIONS; loop_var++)
        {
            for(size_t i = 0; i < ARRAY_SIZE; i++)
            {
                inTemp = in[i];
                outTemp = inTemp * a0 + z1;
                z1 = inTemp * a1 + z2 - b1 * outTemp;
                z2 = inTemp * a2 - b2 * outTemp;
                out[i] = outTemp;
            }
        }
        duration_ms = absolute_time_diff_us(start_time, get_absolute_time()) / 1000;
        printf("[Core #0] Finished iteration #%llu.\n", counter);
        printf("[Core #0] Finished %u samples in %zu milliseconds.\n", ARRAY_SIZE * ITERATIONS, duration_ms);
        printf("[Core #0] %zu kiloSamples per second.\n", ARRAY_SIZE * ITERATIONS / duration_ms);
        counter = counter + 1;
    }
}