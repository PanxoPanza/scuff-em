/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * libhrutil.cc -- library of general-purpose utility routines 
 *
 * homer reid   -- 10/2008
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>

#if defined(_WIN32)
#  include <windows.h>
#  include <process.h>
#else
#  include <sys/resource.h>
#  include <sys/times.h>
#  include <fenv.h>
#endif

#ifdef HAVE_PTHREAD
#  include <pthread.h>
#endif 

#include <math.h>

#include "libhrutil.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#ifdef HAVE_EXECINFO_H
   #include <execinfo.h>
#endif

#ifdef USE_OPENMP
#  include <omp.h>
#endif

#define MAXSTR 1000

/***************************************************************/
/* Timing functions ********************************************/
/***************************************************************/
static double tictoc=0.0;
static unsigned long tictocMem=0;
double Secs()
{ 
  struct timeval tv;
  gettimeofday(&tv, 0);
  return (double)(tv.tv_sec) + 1.0e-6*((double)(tv.tv_usec));
}
void Tic(bool MeasureBytesAllocated)
 { tictoc=Secs(); 
   if (MeasureBytesAllocated) tictocMem=GetMemoryUsage();
 }
double Toc(unsigned long *BytesAllocated)
 { if (BytesAllocated)
    *BytesAllocated=GetMemoryUsage() - tictocMem;
   return Secs()-tictoc; 
 }

/***************************************************************/
/* String functions (not thread-safe) **************************/
/***************************************************************/
/* given "/home/homer/work/MyFile.dat", return "MyFile.dat" */
char *RemoveDirectories(char *s)
{ static char buffer[MAXSTR];
  char *p;

  p=strrchr(s,'/');

  if (!p)
   strncpy(buffer,s,MAXSTR); 
  else
   strncpy(buffer,p+1,MAXSTR);

  return buffer;
}

/* given "/home/homer/work/MyFile.dat", return "/home/homer/work/MyFile" */
char *RemoveExtension(char *s)
{ static char buffer[MAXSTR];
  char *p;

  strncpy(buffer,s,MAXSTR);
  if ( (p=strrchr(buffer,'.')) )
   *p=0;

  return buffer;

}

/* given "/home/homer/work/MyFile.dat", return "dat" */
char *GetFileExtension(char *s)
{
  static char buffer[MAXSTR];
  char *p;

  if ( (p=strrchr(s,'.')) )
   strncpy(buffer,p+1,999);
  else
   buffer[0]=0;

  return buffer;
}

/* given "/home/homer/work/MyFile.dat", return "MyFile" */
char *GetFileBase(const char *s)
{  
#define NUMBUFS 4
  static int bufhead=0;
  static char buffers[NUMBUFS][MAXSTR];
  bufhead = (bufhead+1)%NUMBUFS;
  char *buffer = buffers[bufhead];
  strncpy(buffer, s, MAXSTR);
  char *p = strrchr(buffer,'/');
  char *ret = p ? p+1 : buffer;
  if ( (p=strrchr(ret,'.')) ) *p=0;         // remove file extension
  return ret;
}

/*
 * public domain strtok_r() by Charlie Gordon
 *
 *   from comp.lang.c  9/14/2007
 *
 *      http://groups.google.com/group/comp.lang.c/msg/2ab1ecbb86646684
 *
 *     (Declaration that it's public domain):
 *      http://groups.google.com/group/comp.lang.c/msg/7c7b39328fefab9c
 * 
*/
char* strtok_r_PublicDomain(
    char *str, 
    const char *delim, 
    char **nextp)
{
    char *ret;

    if (str == NULL)
    {
        str = *nextp;
    }

    str += strspn(str, delim);

    if (*str == '\0')
    {
        return NULL;
    }

    ret = str;

    str += strcspn(str, delim);

    if (*str)
    {
        *str++ = '\0';
    }

    *nextp = str;

    return ret;
}

// MS uses "strtok_s" but this is not supported by MinGW apparently,
// so just punt and use the non re-entrant version
// #define strtok_r(s,d,p) strtok(s,d)
// update 20150802 no, this seems to cause core dumps on
// window. instead use the public-domain strtok implementation
// included above.
#if defined(_WIN32)
  #define strtok_r(s,d,p) strtok_r_PublicDomain(s,d,p)
#endif

// portable wrapper around StrCaseCmp
int StrCaseCmp(const char *s1, const char *s2)
{
#if defined(_WIN32)
  return _stricmp(s1,s2);
#else
  return strcasecmp(s1,s2);
#endif
}

/* given "Now is the time for all", return Tokens[0]="Now", */
/* Tokens[1]="is", etc.                                     */
/* Note that s is modified on return.                       */
sVec Tokenize(char *s, const char *Separators)
{ 
  sVec Tokens;
  char *saveptr=0, *Token=strtok_r(s, Separators, &saveptr);
  while (Token)
   { Tokens.push_back(Token);
     Token=strtok_r(0, Separators, &saveptr);
   }
  return Tokens;
}

// old calling convention retained for backward compatibility
int Tokenize(char *s, char **Tokens, int MaxTokens, const char *Separators)
{ sVec TokenVec=Tokenize(s,Separators);
  for(int nt=0; nt<((int)(TokenVec.size())); nt++)
   { Tokens[nt]=TokenVec[nt];
     if (nt==MaxTokens-1) break;
   }
  return TokenVec.size();
}

/***************************************************************/
/* Error-checked vsnprintf *************************************/
/***************************************************************/
int vsnprintfEC(char *str, size_t size, const char *format, va_list ap)
{
  int n = vsnprintf(str, size, format, ap);
  if (unsigned(n) >= size)
    ErrExit("string of length %d exceeded maximum length %lu\n",
	    n, (unsigned long) size);
  return n;
}

/***************************************************************/
/* extension of fopen that accepts an optional colon-separated */
/* search path.                                                */
/* If Path is null this reduces to just the usual fopen.       */
/* If WhichDir is non-null, then on return *WhichDir points to */
/* a static buffer containing the directory in which the file  */
/* was found.                                                  */
/***************************************************************/
FILE *fopenPath(const char *Path, const char *FileName,
                const char *Mode, char **WhichDir)
{ 
  static char DirFound[MAXSTR];
  if (WhichDir) *WhichDir=DirFound;

  if ( FILE *f=fopen(FileName,Mode) )
   { if (strchr(FileName,'/'))
      { snprintf(DirFound,MAXSTR,FileName);
        *strrchr(DirFound,'/') = 0;
      }
     else 
      sprintf(DirFound,".");
     return f;
   }

  const char *p=Path;
  if (!p) return 0;
  while(*p)
   {
     int n=0;
     char FullFileName[MAXSTR];
     while( *p && *p!=':' && n<MAXSTR )
      FullFileName[n++] = *p++;
     strncpy(DirFound,FullFileName,n);
     DirFound[n]=0;
     FullFileName[n++]='/';
     strncpy(FullFileName+n,FileName,MAXSTR-n);
     if ( FILE *f=fopen(FullFileName,Mode) )
      return f;
     if (*p==':') p++;
   }
  return 0;
}

/***************************************************************/
/* Vararg versions of common functions *************************/
/***************************************************************/
FILE *vfopen(const char *format, const char *mode, ...)
{
  COMPLETE_VARARGS2(format, mode, buffer);
  return fopen(buffer,mode);
}

int vmkdir(const char *format, ...)
{
  COMPLETE_VARARGS(format,buffer);
#if defined(_WIN32)
  return _mkdir(buffer);
#else
  return mkdir(buffer,0755);
#endif
}

int vsystem(const char *format, ...)
{
  COMPLETE_VARARGS(format, buffer);
  return system(buffer);
}

// note: unlike the usual strncat, the parameter buflen here is the 
// total number of bytes that will be written to s, *including* what
// is already present in s, and also including the trailing 0.
// Example: 
//  s="hello" + trailing 0
//  buflen=10;
//  format + ...  evaluates to "world"
// --> on return, s="helloworl" + trailing 0
//
void vstrncat(char *s, size_t n, const char *format, ...)
{
  COMPLETE_VARARGS(format,buffer);
  size_t OldLength = strlen(s) + 1;
  if (OldLength>=n) return;
  strncat(s, buffer, n - OldLength);
}

char *vstrappend(char *s, const char *format, ...)
{ COMPLETE_VARARGS(format,buffer);
  if (s==0) 
   s=strdupEC(buffer);
  else
   { int NS=strlen(s), NB=strlen(buffer);
     s = (char *)reallocEC(s, NS+NB+1);
     strcpy(s + NS, buffer);
   }
  return s;
}

char *vstrdup(const char *format, ...)
{
  COMPLETE_VARARGS(format,buffer);
  return strdupEC(buffer);
}

void vsetenv(const char *VariableName, const char *format, ...)
{
  COMPLETE_VARARGS(format,buffer);
#if defined(_WIN32)
  SetEnvironmentVariable(VariableName, buffer);
#else
  setenv(VariableName, buffer, 1);
#endif
}

void ErrExit(const char *format, ...)
{
  va_list ap; 
  char buffer[2*MAXSTR];
  va_start(ap,format);
  vsnprintfEC(buffer,2*MAXSTR,format,ap);
  va_end(ap);
  fprintf(stderr,"error: %s (aborting)\n",buffer);
  Log("error: %s (aborting)",buffer);
  exit(1);
}

void Warn(const char *format, ...)
{
  COMPLETE_VARARGS(format, buffer);
  fprintf(stderr,"**warning: %s \n",buffer);
  Log("warning: %s \n",buffer);
}

/***************************************************************/
/* Simple general-purpose status logging (not thread-safe) *****/
/***************************************************************/
static char *LogFileName=0;
static int LogToConsole=0;

void SetConsoleLogging()
 { LogToConsole=1; }

void SetLogFileName(const char *format, ...)
{ COMPLETE_VARARGS(format, buffer);
  LogToConsole=0;
  if (LogFileName) free(LogFileName);
  LogFileName=strdupEC(buffer);
}

void DoLog(bool ByThread, bool WriteCR, const char *Message)
{
  FILE *f=0;
  if (LogToConsole) 
   f=stdout;
  else if (LogFileName==0) 
   return;
  else
   f = (ByThread ? vfopen("%s.%i","a",LogFileName,GetThreadNum()) : fopen(LogFileName,"a")); 
  if (f==0) return;
  fprintf(f,"%s: %s%s",GetTimeString(),Message,WriteCR ? "\n" : "");
  if (!LogToConsole) fclose(f);
}

void Log(const char *format, ...)
{ COMPLETE_VARARGS(format,buffer);
  DoLog(false,true,buffer);
}

void LogC(const char *format, ...)
{ COMPLETE_VARARGS(format,buffer);
  DoLog(false,false,buffer);
}

void LogThread(const char *format, ...)
{ COMPLETE_VARARGS(format,buffer);
  DoLog(true,true,buffer);
}

void InitializeLog(char *argv0, const char *Path)
{ char *CodeName=GetFileBase(argv0);
  if (LogFileName==0)
   Path ? SetLogFileName("%s/%s.log",Path,CodeName)
        : SetLogFileName("%s.log",CodeName);
  Log("%s running on %s:%d (%s)",CodeName, GetHostName(), getpid(), GetTimeString());
}

// 20120225 thread-safe logging
#ifdef USE_PTHREAD 
static pthread_mutex_t LogMutex = PTHREAD_MUTEX_INITIALIZER;
#endif
void MutexLog(const char *format, ... )
{ COMPLETE_VARARGS(format,buffer);

#ifdef USE_PTHREAD 
  pthread_mutex_lock(&LogMutex);
#endif
  Log(buffer);
#ifdef USE_PTHREAD 
  pthread_mutex_unlock(&LogMutex);
#endif
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void LogPercent(int n, int N, int Gradations)
{ 
  for(int PerCent=0; PerCent<Gradations; PerCent++)
   if ( n == (N*PerCent/Gradations) )
    MutexLog("  %i %% (%i/%i)...",PerCent*100/Gradations,n,N);
}

/***************************************************************/
/* figure out how many processors (cores) are on this machine. */
/***************************************************************/
int GetNumProcs()
{ 
#if defined(USE_OPENMP)
  return omp_get_num_procs();
#elif defined(_WIN32)
  SYSTEM_INFO sysinfo;
  GetSystemInfo(&sysinfo);
  return sysinfo.dwNumberOfProcessors;
#else
  return sysconf(_SC_NPROCESSORS_ONLN);
#endif
  
#if 0 // old non-portable method
  int n,NumProcs;
  FILE *f;

  f=popen("cat /proc/cpuinfo | grep processor | wc -l","r");
  if (!f)
   return 1;

  n=fscanf(f,"%i",&NumProcs);
  if ( n!=1 || NumProcs<=0 )
   return 1;

  return NumProcs;
#endif
} 


/***************************************************************/
/* figure out default number of threads to use.                */
/*                                                             */
/* The static library variable NumThreads and the SetNumThreads*/
/* routine provide a mechanism for API programmers using       */
/* pthreads to specify a number of threads to use. By default  */
/* this variable is unset (set to -1), in which case           */
/* GetNumThreads() will return the number of processor cores.  */
/***************************************************************/
static int NumThreads=-1;
void SetNumThreads(int pNumThreads)
{ NumThreads=pNumThreads; }

int GetNumThreads()
{
#if defined(USE_OPENMP)
  if (NumThreads!=-1) return NumThreads; // F. Ramirez (2021)
  return omp_get_max_threads();
#elif defined(USE_PTHREAD)
  if (NumThreads!=-1) return NumThreads;
  return GetNumProcs();
#else
  return 1;
#endif
}

int GetThreadNum()
{
#if defined(USE_OPENMP)
  return omp_get_thread_num();
#else
  return 0;
#endif
}

/***************************************************************/
/* return 1 if we are running on a core i7 machine (because in */
/* case we will want to set the number of threads used by the  */
/* multithreaded BLAS library equal to half the actual number  */
/* of cores detected by GetNumProcs() because hyperthreading   */
/* screws up the BLAS cache optimizations)                     */
/***************************************************************/
int DetectCorei7()
{ 
  int GrepReturnValue;

  GrepReturnValue=system("cat /proc/cpuinfo | grep -iq i7");
  if(GrepReturnValue==0)
   return 1;
  else
   return 0;
} 

/***************************************************************/
/* nonportable linux memory usage function since getrusage() is*/
/* not useful on linux systems                                 */
/* if MemoryUsage is non-null, it should point to a buffer     */
/* with room for at least 7 unsigned longs, which are filled   */
/* in on return with the various memory usage statistics       */
/***************************************************************/
unsigned long GetMemoryUsage(unsigned long *MemoryUsage)
{
#if defined(_WIN32)
  for (int i=0; i<7; ++i) MemoryUsage[i] = 0;
#else
  FILE *f;
  if ( (f=fopen("/proc/self/statm","r")) )
   { 
     char Line[MAXSTR];
     char *result=fgets(Line,MAXSTR,f);
     fclose(f);
     if (result==0) return 0;

     unsigned long Mem[7];
     sscanf(Line,"%lu %lu %lu %lu %lu %lu %lu",
                  Mem+0, Mem+1, Mem+2, Mem+3, Mem+4, Mem+5, Mem+6);

     long sz = sysconf(_SC_PAGESIZE);
     if (MemoryUsage)
      for(int n=0; n<7; n++)
       MemoryUsage[n] = sz*Mem[n];

     return sz*Mem[0];

   };
#endif
  return 0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void SetCPUAffinity(int WhichCore)
{ 
#if defined(_GNU_SOURCE) && defined(USE_PTHREAD)
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(WhichCore,&cpuset); 
  pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
#else
  (void) WhichCore; // unused
#endif
}

void GetCPUAffinity()
{
#if defined(_WIN32)
#elif 0 // nonportable
  cpu_set_t cpuset;
  sched_getaffinity(0, sizeof(cpu_set_t), &cpuset);
  printf("** CPUAffinity: ");
  for(int i=0; i<8; i++)
   printf("%i", CPU_ISSET(i,&cpuset) ? 1 : 0);
  printf("\n");
#endif
}

/***************************************************************/
/* wait for keypress *******************************************/
/***************************************************************/
void KeyPause()
{ char buffer[100];
  printf("press any key to continue...");
  fgets(buffer,100,stdin);
}

/***************************************************************/
/* relative difference *****************************************/
/***************************************************************/
double RD(double x, double y)
{ 
  if ( x==0.0 && y==0.0 )
   return 0.0;
  return 2.0*fabs(x-y)/(fabs(x)+fabs(y));
}

double RD(cdouble x, cdouble y)
{ 
  double xMag2=norm(x);
  double yMag2=norm(y);
  cdouble z=x-y;
  double zMag2=norm(z);

  if ( z == cdouble(0.0,0.0) )
   return 0.0;
  return 2.0*sqrt(zMag2 / (xMag2+yMag2) );
}

/***************************************************************/
/* error-checked malloc and realloc.                           */
/* note: this malloc() zeros out the buffer before returning,  */
/* i.e. it behaves like calloc().                              */
/***************************************************************/
void *mallocEC(size_t size)
{
  if (size==0)
   return 0;
  void *v=malloc(size);
  if (size && v==0)
   ErrExit("out of memory");
  memset(v,0,size); 
  return v;
}

void *reallocEC(void *v, size_t size)
{
  v=realloc(v, size);
  if (size && v==0)
   ErrExit("out of memory");
  return v;
}

/***************************************************************/
/* error-checked strdup ****************************************/
/***************************************************************/
char *strdupEC(const char *s)
{
  char *dup = (char *) mallocEC(strlen(s) + 1);
  strcpy(dup, s);
  return dup;
}

/***************************************************************/
/* strdup for general pointers *********************************/
/***************************************************************/
void *memdup(void *v, size_t size)
{ if (v==0 || size==0)
   return 0;
  void *vv=mallocEC(size);
  if (vv) memcpy(vv,v,size);
  return vv;
}

/***************************************************************/
/* some complex-number functions *******************************/
/***************************************************************/
cdouble expi(double x)   { return cdouble(cos(x),sin(x)); }

/***************************************************************/
/* complex sqrt that always returns a value with imag part >=0 */
/***************************************************************/
cdouble csqrt2(cdouble z)
{ 
  cdouble rtz=sqrt(z);
  if ( imag(rtz) < 0.0 )
   return -rtz;
  else
   return rtz;
}

/***************************************************************/
/* try to convert string to complex number.                    */
/*   examples: 3+4i, -545.234, -545.234I, 18e3-19e5i           */
/* returns 0 on success, 1 on error.                           */
/***************************************************************/
int S2CD(const char *str, cdouble *pZ)
{ 
  double X, Y;
  char c1, c2;

  int nConv=sscanf(str,"%le%c%le%c",&X,&c1,&Y,&c2);

  if (nConv==1)
   { *pZ=cdouble(X,0.0); return 0; };
  if ( (nConv==2) && (tolower(c1)=='i' || tolower(c1)=='j') )
   { *pZ=cdouble(0.0,X); return 0; };
  if (nConv==2) 
   { *pZ=cdouble(X,0.0); return 0; };
  if (nConv==4 && c1=='+' && (tolower(c2)=='i' || tolower(c1)=='j') )
   { *pZ=cdouble(X,Y); return 0; };
  if (nConv==4 && c1=='-' && (tolower(c2)=='i' || tolower(c1)=='j') )
   { *pZ=cdouble(X,-Y); return 0; };
  return 1;

} 

/***************************************************************/
/* convert complex number to a string                          */
/***************************************************************/
static char DefaultCD2SFormat[100];
static int DefaultCD2SFormatInitialized=0;

char *CD2S(cdouble z, const char *format)  
{ 
  #define NUMSTRINGS 100
  #define STRLEN 50
  static char Strings[NUMSTRINGS][50];
  static int BufPtr=0;

  BufPtr = (BufPtr+1) % NUMSTRINGS;
  char *zStr=Strings[BufPtr];

  if (!format)
   { if ( imag(z)==0.0 )
      snprintf(zStr,MAXSTR,"%g",real(z));
     else if ( real(z)==0.0 )
      snprintf(zStr,MAXSTR,"%gi",imag(z));
     else 
      snprintf(zStr,MAXSTR,"%g+%gi",real(z),imag(z));
     return zStr; 
   };

  int nc, nConv;
  for(nConv=0, nc=0; format[nc]; nc++)
   if (format[nc]=='%') nConv++;

  if (nConv==2)
   snprintf(zStr,STRLEN,format,real(z),imag(z));
  else if (nConv==1)
   snprintf(zStr,STRLEN,format,real(z));
  else 
   strncpy(zStr,format,STRLEN);

  return zStr;
}

char *CD2S(cdouble z)
{ 
  if (DefaultCD2SFormatInitialized==0)
   { strcpy(DefaultCD2SFormat,"(%+9.5e,%+9.5e)");
     DefaultCD2SFormatInitialized=1;
   };
  return CD2S(z,DefaultCD2SFormat);
}

void SetDefaultCD2SFormat(const char *format)
{ 
  strncpy(DefaultCD2SFormat,format,99);
  DefaultCD2SFormatInitialized=1;
}

char *z2s(cdouble z)
{ 
  return CD2S(z,0);
}

void z2s(cdouble z, char *zStr)
 { strcpy(zStr, z2s(z)); }

/***************************************************************/
/* given a string like 'MyFile.dat', try to create a new file  */
/* called 'MyFile.dat.' if that doesn't work (because there is */
/* already a file called 'MyFile.dat', then try to create a    */
/* new file called 'MyFile.dat.1'. if that doesn't work, try   */
/* to create a new file called 'MyFile.dat.2', and so on.      */
/*                                                             */
/* on success, return the new file, opened for writing.        */
/*                                                             */
/* if PrintMessage is nonzero, then we print a console message */
/* in cases in which we had to append an integer to base to    */
/* obtain a unique filename.                                   */
/*                                                             */
/* if FileName is nonzero, then FileName must point to a buffer*/
/* with space for at least MAXSTR characters, which is filled  */
/* in with the name of the file that was eventually created.   */
/***************************************************************/
FILE *CreateUniqueFile(const char *Base, int ConsoleMessage, char *FileName)
 {
   FILE *f;
   char NewFileName[MAXSTR];
   int N;

   for(N=0; N<MAXSTR; N++)
    { 
      if (N==0)
       snprintf(NewFileName,MAXSTR-10,"%s",Base);
      else
       snprintf(NewFileName,MAXSTR-10,"%s.%i",Base,N);

      f=fopen(NewFileName,"r");
      if (f==0) break;
      fclose(f);

    };
   if (f!=0)
    return 0;

   f=fopen(NewFileName,"w");

   if ( ConsoleMessage && N>0 )
    fprintf(stderr,"\n ** WARNING: file %s exists (writing to file %s instead)\n",
            Base,NewFileName);
                  
   if (FileName)
    strncpy(FileName,NewFileName,MAXSTR-1);

   return f;

 } 

// return a string form of the current time. If no buffer is supplied,
// the return value is a pointer to a static buffer, so not thread-safe.
char *GetTimeString(char *s)
{
  static char TimeStr[30];
  if (s==0) s=TimeStr;
  time_t MyTime=time(0);
  strftime(s,30,"%D::%T",localtime(&MyTime));
  return s;
}

char *GetHostName()
{ 
  static char HostName[MAXSTR];
#if defined(_WIN32)
  // gethostname exists on Win32 but causes link failure with MinGW
  // ... punt, since scuff doesn't actually use this function.
  strcpy(HostName, "Windows");
#else
  gethostname(HostName,MAXSTR);
#endif
  return HostName;
}
  

FILE *CreateUniqueFile(const char *Base, int ConsoleMessage) 
 { return CreateUniqueFile(Base,ConsoleMessage,0); }

FILE *CreateUniqueFile(const char *Base) 
 { return CreateUniqueFile(Base,0,0); }

/***************************************************************/
/***************************************************************/
/***************************************************************/

#if !defined(_WIN32)
static char CodeMarker[MAXSTR];
void SetCodeMarker(const char *Marker)
{ strncpy(CodeMarker, Marker, MAXSTR); }

static void SignalHandler(int WhichSignal)
{
  FILE *f = (LogFileName ? fopen(LogFileName,"a") : stderr );
  fprintf(f,"\n\n**<><> received signal %i ",WhichSignal);
  if (CodeMarker[0]!=0)
   fprintf(f,"while %s",CodeMarker);
  fprintf(f,"\n\n");

#ifdef HAVE_BACKTRACE
  void *buffer[20];
  int nptrs=backtrace(buffer, 20);
  backtrace_symbols_fd(buffer, nptrs, fileno(f));
#endif

  fclose(f);
  exit(1);
}
#endif // !defined(_WIN32)

void InstallHRSignalHandler()
{ 
#if defined(_WIN32) || defined(__APPLE__)
#else

  // this allows users to set the environment variable
  //  SCUFF_ABORT_ON_FPE=1
  // to force the code to abort on any floating-point error
  char *s=getenv("SCUFF_ABORT_ON_FPE");
  if (s && s[0]=='1')
   feenableexcept(FE_INVALID | FE_OVERFLOW);

  struct sigaction sa;

  sa.sa_handler=SignalHandler;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags=0;
  //sigaction(SIGKILL, &sa, 0); //can't mask this one 
  sigaction(SIGSEGV, &sa, 0);
  sigaction(SIGPIPE, &sa, 0);
  sigaction(SIGALRM, &sa, 0);
  sigaction(SIGFPE,  &sa, 0);
  sigaction(SIGTRAP, &sa, 0);
  sigaction(SIGABRT, &sa, 0);
  //sigaction(SIGCHLD, &sa, 0);

  CodeMarker[0]=0;
#endif
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool IsFinite(double d)
{ return !(ISNAN(d));
#if 0
// 20120622 I don't understand why fpclassify() seems not to 
//          be portable; isn't it a POSIX function? I am 
//          removing this for now. 
  // int fpc = fpclassify(d);
  int fpc = __fpclassify(d);
  if (fpc==FP_NORMAL || fpc==FP_ZERO)
   return 1;
  return 0;
#endif
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool IsFinite(cdouble z)
 { return IsFinite(real(z)) && IsFinite(imag(z)); }
 
/***************************************************************/
/***************************************************************/
/***************************************************************/
void PlusEqualsVec(double *V1, double Alpha, double *V2, int N)
 { for(int n=0; n<N; n++) V1[n]+=Alpha*V2[n]; }
void PlusEqualsVec(cdouble *V1, cdouble Alpha, cdouble *V2, int N)
 { for(int n=0; n<N; n++) V1[n]+=Alpha*V2[n]; }
void PlusEqualsVec(cdouble *V1, cdouble Alpha, double *V2, int N)
 { for(int n=0; n<N; n++) V1[n]+=Alpha*V2[n]; }
void PlusEqualsVec(double *V1, double *V2, int N)
 { for(int n=0; n<N; n++) V1[n]+=V2[n]; }
void PlusEqualsVec(cdouble *V1, cdouble *V2, int N)
 { for(int n=0; n<N; n++) V1[n]+=V2[n]; }
void ScaleVec(double *V1, double Alpha, int N)
 { for(int n=0; n<N; n++) V1[n]*=Alpha; }
void ScaleVec(cdouble *V1, cdouble Alpha, int N)
 { for(int n=0; n<N; n++) V1[n]*=Alpha; }

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool CheckEnv(const char *Name, const char *fmt, void *Destination, bool LogSuccess)
{ char *Value=getenv(Name);
  bool Status=(Value && 1==sscanf(Value,fmt,Destination));
  if (Status && LogSuccess) Log("Set %s to %s",Name,Value);
  return Status;
}
bool CheckEnv(const char *Name, double *Destination, bool LogSuccess)
{ return CheckEnv(Name, "%le", (void *)Destination, LogSuccess); }
bool CheckEnv(const char *Name, int *Destination, bool LogSuccess)
{ return CheckEnv(Name, "%i", (void *)Destination, LogSuccess); }

bool CheckEnv(const char *Name, char **Destination, bool LogSuccess)
{ char *Value = getenv(Name);
  if (Value)
   { *Destination = strdup(Value);
     if (LogSuccess) Log("Set %s to %s",Name,Value);
     return true;
   }
  return false;
}

bool CheckEnv(const char *Name, bool LogSuccess)
{ char *s=getenv(Name);
  bool Status = (s && s[0]=='1');
  if (Status && LogSuccess) Log("Set %s = true",Name);
  return Status;
}

void AppendEnv(const char *Name, const char *Separator, const char *format, ...)
{ 
  COMPLETE_VARARGS(format, buffer);
  char *Old = getenv(Name);
  if (Old) { // check if Old already contains buffer
    const char *p = strstr(Old, buffer);
    if (p) {
      size_t seplen = strlen(Separator);
      size_t buflen = strlen(buffer);
      if ((p == Old || (p >= Old + seplen && !strncmp(p-seplen, Separator, seplen)))
          && (*(p+buflen)==0 || !strcmp(p+buflen, Separator)))
        return; // Old contains buffer with Separator before & after
    }
  }
  char *New = vstrdup("%s%s%s",Old ? Old:"",Old&&Separator ? Separator:"",buffer);
  setenv(Name,New,1);
  free(New);
  return;
}
