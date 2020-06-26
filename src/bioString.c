#include <string.h>
#include <ctype.h>

char* rstrip(char* s)
{
    char* p = s + strlen(s);
    while (p > s && isspace((unsigned char) (*--p)))
        *p = '\0';
    return s;
}

char* stripComment(char* s)
{
    char* p = s;
    while (*p != '!' && *p != '\0')
    {
        p++;
    }
    *p = '\0';
    return s;
}
