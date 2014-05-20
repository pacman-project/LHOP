
#pragma once
#ifndef _PLATFORM_
#define _PLATFORM_

#if defined WIN32 | defined WIN64
#define FILE_SEPARATOR '\\'
#define TEMPORARY_DIRECTORY "C:\\Temp\\"
#else
#define FILE_SEPARATOR '/'
#define TEMPORARY_DIRECTORY "/tmp/"
#endif


#endif


