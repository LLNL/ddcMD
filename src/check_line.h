#ifndef CHECK_LINE_H
#define CHECK_LINE_H
#include "pio.h"
#include "gid.h"

char *checkRecord(int *error, unsigned recordLength, char *record, PFILE *pfile);
void readline_init(PFILE *pfile);
int readline(char **tail_ptr, PFILE *pfile);
int parsehead(char *line, char *labelFmt, gid_type* label, char **class,
	      char **tail, int checksum_ok);
int parseatom(char *line, char **tail_ptr, char **name_ptr,
	      char **group_ptr, double r[3], double v[3], int checksum_ok);

#endif


/* Local Variables: */
/* tab-width: 3 */
/* End: */
