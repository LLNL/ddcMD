#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "sCatPrintf.h"

// requirement: *line has to be NULL or a pointer to address of a string

int sCatPrintf(char **line, char *format, ...)
{
    va_list args;
    va_start(args, format);
    char *new;
    int sizeNew = vasprintf(&new, format, args);
    va_end(args);
    if (sizeNew < 0) return sizeNew;
    if (*line == NULL) *line = calloc(0, 1); //if *line is NULL generate a NULL string
    int sizeLine = strlen(*line);
    sizeLine += sizeNew + 1;
    *line = realloc(*line, sizeLine);
    if (line == NULL) return -2;
    strcat(*line, new);
    free(new);
    return sizeLine - 1; // size of line excluding the terminating null byte ('\0')
}
