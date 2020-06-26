#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#define xstr(x) str(x)
#define str(x...) #x
#define CompileDate  __DATE__  
#define CompileTime  __TIME__  
int getRank(int);
static char* devlpmtSrc = NULL;
static char* srcpath = NULL;
static char* cflags = NULL;
static char* ldflags = NULL;
static char* arch = NULL;
static char* target = NULL;
static char* svnversion = NULL;
static char* executable = NULL;
static char compiletime[64];
static char version[128];

static void  version_print(FILE *file, int argc, char *argv[])
{

    char hostname[256],cwd[256];
    gethostname(hostname, sizeof(hostname));
    getcwd(cwd, 255);
    printf("ddcMD run info: MPIRANK=%d HOSTNAME=%s PID=%d CWD=%s\n", getRank(0),hostname,getpid(), cwd); fflush(stdout);
#	  if defined(BANNER)
         fprintf(file,BANNER, svnversion);
#	  endif
	  fprintf(file,"cmdLine    ="); 
	  for (int i=0;i<argc;i++) fprintf(file,"%s ",argv[i]); 
      fprintf(file,"\n"); 
      fprintf(file,"version    =%s\n", version);
      fprintf(file,"git_version=%s\n", svnversion);
      fprintf(file,"srcpath    =%s\n", srcpath);
      fprintf(file,"cflags     =%s\n", cflags);
      fprintf(file,"ldflags    =%s\n", ldflags);
      fprintf(file,"compiletime=%s\n", compiletime);
      fprintf(file,"target     =%s\n", target);
      fprintf(file,"executable =%s\n", executable );
      fprintf(file,"devlpmtSrc =\n %s",devlpmtSrc); 
}
void version_init(int argc, char *argv[])
{
   srcpath    = strdup(xstr(PATH));
   cflags     = strdup(xstr(CFLAGS));
   ldflags    = strdup(xstr(LDFLAGS));
   arch       = strdup(xstr(ARCH));
   target     = strdup(xstr(TARGET));
   svnversion = strdup(xstr(GITVERSION));
   executable = strdup(argv[0]);
   sprintf(compiletime,"%s %s",CompileDate,CompileTime);
   sprintf(version,"%s r%s (Compiled on %s)",target, svnversion, compiletime);
   devlpmtSrc = strdup(xstr(DEVLPMT_SRC));
   for (unsigned i=0;i<strlen(devlpmtSrc)+1;i++) if (devlpmtSrc[i] == ';') devlpmtSrc[i]='\n'; 
   if (getRank(0)==0)
   {
		version_print(stdout, argc, argv); 
		FILE *file = fopen("ddcMD.header","w"); 
		version_print(file, argc, argv); 
		fclose(file); 
   }
}
char *version_get(char *string)
{
   if (!strcmp(string, "srcpath"))     return srcpath; 
   if (!strcmp(string, "cflags"))      return cflags; 
   if (!strcmp(string, "ldflags"))     return ldflags; 
   if (!strcmp(string, "arch"))        return arch; 
   if (!strcmp(string, "target"))      return target; 
   if (!strcmp(string, "svnversion"))  return svnversion; 
   if (!strcmp(string, "executable"))  return executable; 
   if (!strcmp(string, "compiletime")) return compiletime; 
   if (!strcmp(string, "version"))     return version; 
   assert(1==0);
   return NULL; 
}


/* Local Variables: */
/* tab-width: 3 */
/* End: */
