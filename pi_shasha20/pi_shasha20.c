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
#define         ITERATIONS      256

#define     U32_LEFT_ROTATE(input, dist)        (((input) << dist) | ((input) >> (32u - dist)))
#define     U32_RIGHT_ROTATE(input, dist)       (((input) >> dist) | ((input) << (32u - dist)))
#define     LOAD_U32_BE(buffer)                 ((uint32_t)(buffer)[3u] | ((uint32_t)(buffer)[2u] << 8u) | ((uint32_t)(buffer)[1u] << 16u) | ((uint32_t)(buffer)[0u] << 24u))   
#define     LOAD_U32_LE(buffer)                 ((uint32_t)(buffer)[0u] | ((uint32_t)(buffer)[1u] << 8u) | ((uint32_t)(buffer)[2u] << 16u) | ((uint32_t)(buffer)[3u] << 24u)) 
#define     STORE_U32_BE(buffer, input)         for(size_t macro_loop_var = 0u; macro_loop_var < 4u; macro_loop_var++)\
                                                {\
                                                    (buffer)[macro_loop_var] = (input) >> (8u * (3u - macro_loop_var));\ 
                                                }             
#define     CHACHA20_QUARTER_ROUND(a, b, c, d)  a += b; \
                                                d = U32_LEFT_ROTATE(d ^ a, 16u); \
                                                c += d; \
                                                b = U32_LEFT_ROTATE(b ^ c, 12u); \
                                                a += b; \
                                                d = U32_LEFT_ROTATE(d ^ a,  8u); \
                                                c += d; \
                                                b = U32_LEFT_ROTATE(b ^ c,  7u)

static const uint32_t SHA256_constants[64u] =
{   0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
    0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u, 0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
    0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
    0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
    0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u, 0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
    0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
    0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
    0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u, 0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u };

static void SHA256_state_init(uint32_t* SHA256_output)
{
    SHA256_output[0u] = 0x6a09e667u;
    SHA256_output[1u] = 0xbb67ae85u;
    SHA256_output[2u] = 0x3c6ef372u;
    SHA256_output[3u] = 0xa54ff53au;
    SHA256_output[4u] = 0x510e527fu;
    SHA256_output[5u] = 0x9b05688cu;
    SHA256_output[6u] = 0x1f83d9abu;
    SHA256_output[7u] = 0x5be0cd19u;
}

static void SHA256_block_processor(uint8_t* SHA256_working_buffer, uint32_t* SHA256_output)
{
    size_t loop_var;
    uint32_t temp_vars[6u];
    uint32_t round_hash_values[8u];

    for(loop_var = 16u; loop_var < 64u; loop_var++)
    {
        temp_vars[2u] = LOAD_U32_BE(SHA256_working_buffer + (4u * (loop_var - 15u)));
        temp_vars[3u] = LOAD_U32_BE(SHA256_working_buffer + (4u * (loop_var - 2u)));
        temp_vars[0u] = U32_RIGHT_ROTATE(temp_vars[2u], 7u) ^ U32_RIGHT_ROTATE(temp_vars[2u], 18u) ^ (temp_vars[2u] >> 3u);
        temp_vars[1u] = U32_RIGHT_ROTATE(temp_vars[3u], 17u) ^ U32_RIGHT_ROTATE(temp_vars[3u], 19u) ^ (temp_vars[3u] >> 10u);
        STORE_U32_BE(SHA256_working_buffer + (4u * loop_var), LOAD_U32_BE(SHA256_working_buffer + (4u * (loop_var - 16u))) + temp_vars[0u] + LOAD_U32_BE(SHA256_working_buffer + (4u * (loop_var - 7u))) + temp_vars[1u]);
    }

    memcpy(round_hash_values, SHA256_output, sizeof(uint32_t) * 8);

    for(loop_var = 0u; loop_var < 64u; loop_var++)
    {
        temp_vars[0u] = U32_RIGHT_ROTATE(round_hash_values[4u], 6u) ^ U32_RIGHT_ROTATE(round_hash_values[4u], 11u) ^ U32_RIGHT_ROTATE(round_hash_values[4u], 25u);
        temp_vars[1u] = (round_hash_values[4u] & round_hash_values[5u]) ^ (~(round_hash_values[4u]) & round_hash_values[6u]);
        temp_vars[2u] = round_hash_values[7u] + temp_vars[0u] + temp_vars[1u] + SHA256_constants[loop_var] + LOAD_U32_BE(SHA256_working_buffer + (4u * loop_var));
        temp_vars[3u] = U32_RIGHT_ROTATE(round_hash_values[0u], 2u) ^ U32_RIGHT_ROTATE(round_hash_values[0u], 13u) ^ U32_RIGHT_ROTATE(round_hash_values[0u], 22u);
        temp_vars[4u] = (round_hash_values[0u] & round_hash_values[1u]) ^ (round_hash_values[0u] & round_hash_values[2u]) ^ (round_hash_values[1u] & round_hash_values[2u]);
        temp_vars[5u] = temp_vars[3u] + temp_vars[4u];

        round_hash_values[7u] = round_hash_values[6u];
        round_hash_values[6u] = round_hash_values[5u];
        round_hash_values[5u] = round_hash_values[4u];
        round_hash_values[4u] = round_hash_values[3u] + temp_vars[2u];
        round_hash_values[3u] = round_hash_values[2u];
        round_hash_values[2u] = round_hash_values[1u];
        round_hash_values[1u] = round_hash_values[0u];
        round_hash_values[0u] = temp_vars[2u] + temp_vars[5u];
    }

    for(loop_var = 0u; loop_var < 8u; loop_var++)
    {
        SHA256_output[loop_var] = SHA256_output[loop_var] + round_hash_values[loop_var];
    }
}

static void SHA256_eof_processor(uint8_t* SHA256_working_buffer, size_t read_bytes, uint32_t input_size, uint32_t* SHA256_output)
{
    if(read_bytes < 56u)
    {
        SHA256_working_buffer[read_bytes++] = 0x80u;
        for(; read_bytes < 60u; read_bytes++)
        {
            SHA256_working_buffer[read_bytes] = 0u;
        }
    }
    else
    {
        SHA256_working_buffer[read_bytes++] = 0x80u;
        for(; read_bytes < 64u; read_bytes++)
        {
            SHA256_working_buffer[read_bytes] = 0u;
        }
        SHA256_block_processor(SHA256_working_buffer, SHA256_output);

        for(read_bytes = 0u; read_bytes < 60u; read_bytes++)
        {
            SHA256_working_buffer[read_bytes] = 0u;
        }
    }
    
    for(; read_bytes < 64u; read_bytes++)
    {
        SHA256_working_buffer[read_bytes] = input_size >> (8u * (63u - read_bytes));
    }
    SHA256_block_processor(SHA256_working_buffer, SHA256_output);
}

static void chacha20_state_block_init(uint32_t* chacha20_state_block, uint32_t* key)
{
    chacha20_state_block[0u] = LOAD_U32_LE("expa");
    chacha20_state_block[1u] = LOAD_U32_LE("nd 3");
    chacha20_state_block[2u] = LOAD_U32_LE("2-by");
    chacha20_state_block[3u] = LOAD_U32_LE("te k");
    memcpy(chacha20_state_block + 4u, key, sizeof(uint32_t) * 8u);
    chacha20_state_block[14u] = LOAD_U32_LE("mist");
    chacha20_state_block[15u] = LOAD_U32_LE("rake");
}

static void gen_chacha20_xor_block(uint32_t* chacha20_state_block, uint32_t* chacha20_xor_block, uint64_t block_counter)
{
    chacha20_state_block[12u] = block_counter & 0xffffffffu;
    chacha20_state_block[13u] = block_counter >> 32u;

    memcpy(chacha20_xor_block, chacha20_state_block, sizeof(uint32_t) * 16u);

    for (size_t loop_var = 0u; loop_var < 10u; loop_var++) /* 20 rounds, 2 rounds per loop. */
    {
        CHACHA20_QUARTER_ROUND(chacha20_xor_block[0u], chacha20_xor_block[4u], chacha20_xor_block[8u], chacha20_xor_block[12u]);    /* column 0 */
        CHACHA20_QUARTER_ROUND(chacha20_xor_block[1u], chacha20_xor_block[5u], chacha20_xor_block[9u], chacha20_xor_block[13u]);    /* column 1 */
        CHACHA20_QUARTER_ROUND(chacha20_xor_block[2u], chacha20_xor_block[6u], chacha20_xor_block[10u], chacha20_xor_block[14u]);   /* column 2 */
        CHACHA20_QUARTER_ROUND(chacha20_xor_block[3u], chacha20_xor_block[7u], chacha20_xor_block[11u], chacha20_xor_block[15u]);   /* column 3 */
        CHACHA20_QUARTER_ROUND(chacha20_xor_block[0u], chacha20_xor_block[5u], chacha20_xor_block[10u], chacha20_xor_block[15u]);   /* diagonal 0 */
        CHACHA20_QUARTER_ROUND(chacha20_xor_block[1u], chacha20_xor_block[6u], chacha20_xor_block[11u], chacha20_xor_block[12u]);   /* diagonal 1 */
        CHACHA20_QUARTER_ROUND(chacha20_xor_block[2u], chacha20_xor_block[7u], chacha20_xor_block[8u], chacha20_xor_block[13u]);    /* diagonal 2 */
        CHACHA20_QUARTER_ROUND(chacha20_xor_block[3u], chacha20_xor_block[4u], chacha20_xor_block[9u], chacha20_xor_block[14u]);    /* diagonal 3 */                                
    }

    for(size_t loop_var = 0u; loop_var < 16u; loop_var++) /* adding original block to scrambled block */
    {
        chacha20_xor_block[loop_var] = chacha20_xor_block[loop_var] + chacha20_state_block[loop_var];
    }
}

void shasha20_processor(uint8_t* buffer, size_t buffer_length, size_t iteration_count, int core_number)
{
    absolute_time_t start_time = get_absolute_time();
    uint32_t SHA256_output[8];
    uint8_t SHA256_working_buffer[256];

    for(size_t i = 0; i < iteration_count; i++)
    {
        SHA256_state_init(SHA256_output);
        for(size_t j = 0; j < (buffer_length / 64); j++)
        {
            memcpy(SHA256_working_buffer, buffer + (64 * j), 64);
            SHA256_block_processor(SHA256_working_buffer, SHA256_output);
        }
        SHA256_eof_processor(SHA256_working_buffer, buffer_length % 64, buffer_length * 8, SHA256_output);
    }

    size_t duration_ms = absolute_time_diff_us(start_time, get_absolute_time()) / 1000;
    printf("[Core #%d] Finished %zu SHA256 iterations in %zu milliseconds.\n", core_number, iteration_count, duration_ms);
    printf("[Core #%d] Speed: %zu kilobytes of SHA256 per second.\n", core_number, iteration_count * buffer_length / duration_ms);

    start_time = get_absolute_time();
    uint32_t chacha20_state_block[16u];
    uint32_t chacha20_xor_block[16u];
    uint64_t block_counter = 0u;
    chacha20_state_block_init(chacha20_state_block, SHA256_output);
    for(size_t i = 0; i < iteration_count; i++)
    {
        for(size_t j = 0; j < buffer_length; j = j + 64u)
        {
            gen_chacha20_xor_block(chacha20_state_block, chacha20_xor_block, block_counter);
            block_counter = block_counter + 1;
        }
    }

    gen_chacha20_xor_block(chacha20_state_block, chacha20_xor_block, block_counter);
    for(size_t i = 0; i < (buffer_length % 64u); i++)
    {
        buffer[((buffer_length / 64u) * 64u) + i] = buffer[((buffer_length / 64u) * 64u) + i] ^ (chacha20_xor_block[i / 4u] >> ((i % 4u) * 8u));
    }
    duration_ms = absolute_time_diff_us(start_time, get_absolute_time()) / 1000;
    printf("[Core #%d] Finished %zu ChaCha20 iterations in %zu milliseconds.\n", core_number, iteration_count, duration_ms);
    printf("[Core #%d] Speed: %zu kilobytes of ChaCha20 per second.\n", core_number, iteration_count * buffer_length / duration_ms);
}

void core1_main(void)
{
    uint8_t buffer[65536];
    uint64_t counter = 1;
    for(size_t i = 0; i < sizeof(buffer); i++)
    {
        buffer[i] = i & 0xFF;
    }
    while(1)
    {
        printf("[Core #1] Beginning run #%llu.\n", counter);
        shasha20_processor(buffer, sizeof(buffer), ITERATIONS, 1);
        counter = counter + 1;
    }   
}

int main(void)
{
    // sleep_ms(1000);
    // vreg_set_voltage(VREG_VOLTAGE_0_90);
    // sleep_ms(1000);
    // if(set_sys_clock_khz(32000, false))
    // {
    //     sleep_ms(1000);
    //     gpio_init(PICO_DEFAULT_LED_PIN);
    //     gpio_set_dir(LED_PIN, GPIO_OUT);
    //     gpio_put(LED_PIN, 1);
    // }

    stdio_usb_init();
    multicore_launch_core1(core1_main);
    uint64_t counter = 1;

    uint8_t buffer[65536];
    for(size_t i = 0; i < sizeof(buffer); i++)
    {
        buffer[i] = i & 0xFF;
    }
    while(1)
    {
        printf("[Core #0] Beginning run #%llu.\n", counter);
        shasha20_processor(buffer, sizeof(buffer), ITERATIONS, 0);
        counter = counter + 1;
   }