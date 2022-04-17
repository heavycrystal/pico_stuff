cmake_minimum_required(VERSION 3.12)

# Pull in SDK (must be before project)
include(pico_sdk_import.cmake)

project(pi_shasha20 C CXX ASM)
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)

# Initialize the SDK
pico_sdk_init()

add_compile_options(-Wall
        -Wno-format
        -Wno-unused-function
        -Wno-maybe-uninitialized
        )

add_executable(pi_shasha20
	pi_shasha20.c
)

pico_define_boot_stage2(slower_boot2 /home/kevin/gen_coding/pico_stuff/pico-sdk/src/rp2_common/boot_stage2/compile_time_choice.S)
target_compile_definitions(slower_boot2 PRIVATE PICO_FLASH_SPI_CLKDIV=4)

pico_set_boot_stage2(pi_shasha20 slower_boot2)

target_link_libraries(pi_shasha20 pico_stdlib pico_multicore)

pico_enable_stdio_usb(pi_shasha20 1)
pico_enable_stdio_uart(pi_shasha20 0)

pico_add_extra_outputs(pi_shasha20)