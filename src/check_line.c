#include "check_line.h"

#include <stdlib.h>
#include <sys/types.h>
#include <stdio.h>
#include <regex.h>
#include <string.h>
#include <ctype.h>
#include "three_algebra.h"
#include "object.h"
#include "pio.h"
#include "crc32.h"
#include "species.h"
#include "match.h"

enum CHKLINE {FIXED=1,CHECKSUM=2,LENGTH_WRONG=4,NO_NEWLINE=8,BAD_CHARACTER=16,NOCHECKSUM=32};

int getRank(int);
int checkpos(char *string, double *r);
int checkvel(char *string, double *v);
int checklabel(char *string);

//static PFILE *pfile;
static char name_pattern[1024],type_pattern[1024],*group_pattern="TEMP|FREE"; 
static regex_t *class_pattern; 
static const int
  n_pat_len = sizeof(name_pattern)/sizeof(*name_pattern),
  t_pat_len = sizeof(type_pattern)/sizeof(*type_pattern);

void readline_init(PFILE *pfile)
{
	char **namelist,**typelist;
	species_get(NULL,NAMELIST,(void *)&namelist);
	int nspecies=0;
	name_pattern[0] = '\0';
	strncat(name_pattern,"^(",n_pat_len - strlen(name_pattern) - 1);
	if (namelist[nspecies]!=NULL)
	  strncat(name_pattern,namelist[nspecies++],n_pat_len - strlen(name_pattern) - 1);
	while (namelist[nspecies] != NULL) 
	{
	  strncat(name_pattern,"|",n_pat_len - strlen(name_pattern) - 1);
	  strncat(name_pattern,namelist[nspecies++],n_pat_len - strlen(name_pattern) - 1);
	}
	strncat(name_pattern,") ",n_pat_len - strlen(name_pattern) - 1);
	species_get(NULL,TYPELIST,(void *)&typelist);
	int ntypes=0; 
	type_pattern[0] = '\0';
	if (typelist[ntypes]!=NULL)
	  strncat(type_pattern,typelist[ntypes++],t_pat_len - strlen(type_pattern) - 1);
	while (typelist[ntypes] != NULL) 
	{
	  strncat(type_pattern,"|",t_pat_len - strlen(type_pattern) - 1);
	  strncat(type_pattern,typelist[ntypes++],t_pat_len - strlen(type_pattern) - 1);
	}
	class_pattern = compilepattern("ATOMS?");
	pfile->recordIndex--; 
}
/*
int readline(char **tail_ptr, PFILE *pfile)
{
   int maxRecordLength=2048; 
   char record[maxRecordLength];
   Pfgets(record, maxRecordLength-1, pfile);
   size_t recordLength = strnlen(record,maxRecordLength);
   pfile->recordIndex++; 
}
*/

// The character '!' is not allowed in records. It is used flag bad characters.   All bad characters when founded are replaced by '!'. 
char *checkRecord(int *error, unsigned recordLength, char *record, PFILE *pfile)
{
   *error=0;
   if (pfile->recordLength != 0 && recordLength != pfile->recordLength)  *error |=LENGTH_WRONG;
   if (record[recordLength-1] != '\n') *error |= NO_NEWLINE;
   for (unsigned i=0;i<recordLength-1;i++)  
   {
      if ( record[i]=='\0' ) 
      {
         record[i]='!';   //  '\0'  is a bad character and is replace \0 with '!'.   This will allow for string processing of record.  
         *error |= BAD_CHARACTER ;
      }
   }
   if (pfile->checksum != CRC32) 
   {
      *error |= NOCHECKSUM; 
      return record; 
   }
   char *end; 
   if (record[8] != ' ')    // record[8] must be ' '.  If not try setting it to ' ' and see if checksum works. 
   { 
      record[8]=' '; 
      *error |= FIXED; 
   }
   unsigned long checksum=strtoul(record,&end,16);
   unsigned checksumnew  = checksum_crc32((unsigned char *)record+8,pfile->recordLength-8);
   record[8]='\0';
   if (record+8!=end)
   {
      *end='!';
      *error |= CHECKSUM; 
      return record+9; 
   }
   if (checksum != checksumnew)
   {
         *error |= CHECKSUM; 
         return record+9; 
   }
   return record+9; 
}
int parsehead(char *line, char *labelFmt, gid_type *label, char **class_field, char **tail_field,int checksum_ok)
{
   unsigned long label_offset,label_ok; 
   unsigned long class_ok;
   char *class,*end_class; 
   char *label_field; 

   label_field=strtok(line," ");
   *class_field=strtok(NULL," ");
   *tail_field=*class_field+strlen(*class_field)+1;
   label_offset=0;
   sscanf(label_field, labelFmt,label,&label_offset);
   if  (checksum_ok) return 0; 
   label_ok=class_ok=0; 
   if ( label_field && !index(label_field,'-')) 
   {
      if (label_offset+label_field-line >= 8 && label_offset <= 10 && *(label_field+label_offset) == '\0')  label_ok=1;
      else
         checklabel(label_field);
   }
   class = findpattern1(*class_field,class_pattern,&end_class); 
   /*   	class = findpattern(*class_field,"ATOMS?",&end_class); */
   if (class) 
   {
      if (class == *class_field) 
      {
         class_ok=1;
         if ( class[5] == 'S' ) class[5]= '\0'; 
         *tail_field= end_class + 1; 
      }
   }

   if (class_ok && label_ok) return 0; 
   return BAD_CHARACTER; 
}
int parseatom(char *line,char **tail_ptr,char **name_ptr, char **group_ptr, double r[3], double v[3],int checksum_ok)
{
   int i, error;
   char *tail,*end,*name,*group,*floats,*decimal;
   error=0;
   tail = line; 
   *name_ptr=name = findpattern(line,name_pattern,&end);
   if (name)
   {
      *(end-1)='\0' ;
      tail = end;
   }
   *group_ptr=group = findpattern(tail,group_pattern,&end);
   if (group)
   {
      group[-1]=*(end)='\0' ;
      tail = end+1; 
   }
   decimal = findpattern(tail,"\\.[0-9]{13}",&end);
   floats = NULL; 
   if (decimal)
   {
      floats = decimal-4;
      while (floats  > tail + 19 ) floats-=19;
      end=floats+18;
      if (!group) {*(floats-1)='\0';}
      for(i=0;i<3;i++)
      {
         floats[18+i*19]='\0';
         error |= checkpos(floats+i*19,r+i);
      }
      for(i=3;i<6;i++)
      {
         floats[18+i*19]='\0';
         error |=checkvel(floats+i*19,v+i);
      }
      tail=floats+6*19;
   }
   if (name && group && floats ) return error; 
   error |= BAD_CHARACTER; 
   return error; 
}
int checklabel(char *string)
{
   char *c,*p[]={" 0123456789","0123456789"}; 
   int ip; 
   ip=0; 
   while(*string != '\0')
   {
      c=index(p[ip],*string);
      if (c) { if (*c != ' ') ip=1; }
      else *string='!'; 
      string++; 
   }
   return 0;
}
int checkpos(char *string,double *r)
{
   int i,ip,error; 
   char *c,*p[]={" -"," -0123456789","0123456789"},*tail;
   error=0; 
   if (string[4]!='.' ){string[4]='.' ; error |= FIXED;}
   *r = strtod(string,&tail);
   if (string+18==tail) return error;
   ip = 0; 
   for (i=0;i<18;i++)
   {
      switch (i)
      {
         case 0:
            c=index(p[ip],string[i]);
            ip=1; 
            if (c) {if (*c!=' ') ip=2; }
            else {string[i]='!' ; error |= BAD_CHARACTER;}
            break; 
         case 1:
            c=index(p[ip],string[i]);
            if (c) {if (*c!=' ') ip=2; }
            else {string[i]='!' ; error |= BAD_CHARACTER;}
            break; 
         case 2:
            c=index(p[ip],string[i]);
            ip=2; 
            if (!c) {string[i]='!' ; error |= BAD_CHARACTER;}
            break; 
         case 3:
            c=index(p[2],string[i]);
            if (!c) {string[i]='!' ; error |= BAD_CHARACTER;}
            break; 
         case 4:
            if (string[i]!='.' ){string[i]='.' ; error |= FIXED;}
            break; 
         case 5:
         case 6:
            c=index(p[2],string[i]);
            if (!c) {string[i]='!' ; error |= BAD_CHARACTER;}
            break; 
         default:
            c=index(p[2],string[i]);
            if (!c) {string[i]='0' ; error |= FIXED;}
            break; 
      }
   }
   if  (error <=FIXED) *r=strtod(string,&tail);
   return error;
}
int checkvel(char *string, double *v)
{
   int i,error; 
   char *c,*p[]={" -","0123456789"}, *tail;
   error=0; 
   if (string[4]!='.' ){string[4]='.' ; error |= FIXED;}
   *v = strtod(string,&tail);
   if (string+18==tail) return error;
   for (i=0;i<18;i++)
   {
      switch (i)
      {
         case 0:
         case 1:
            if (string[i]!=' ' ){string[i]=' ' ; error |= FIXED;}
            break; 
         case 2:
            c=index(p[0],string[i]);
            if (!c) string[i]='!' ;
            break; 
         case 3:
         case 5:
            if (string[i]!='0' ){string[i]='0' ; error |= FIXED;}
            break; 
         case 4:
            if (string[i]!='.' ){string[i]='.' ; error |= FIXED;}
            break; 
         case 6:
         case 7:
            c=index(p[1],string[i]);
            if (!c) string[i]='!' ;
            break; 
         default:
            c=index(p[1],string[i]);
            if (!c) {string[i]='0' ; error |= FIXED;}
            break; 
      }
   }
   if  (error <=FIXED) *v=strtod(string,&tail);
   return error;
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
