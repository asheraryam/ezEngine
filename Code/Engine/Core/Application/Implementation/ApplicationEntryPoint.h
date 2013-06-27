
#pragma once

#if EZ_ENABLED(EZ_PLATFORM_WINDOWS)

  #include <Core/Application/Implementation/Win/ApplicationEntryPoint_win.h>

#elif EZ_ENABLED(EZ_PLATFORM_OSX)

  #include <Core/Application/Implementation/Posix/ApplicationEntryPoint_posix.h>

#else
  #error "Missing definition of platform specific entry point!"
#endif