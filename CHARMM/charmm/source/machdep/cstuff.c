// *********************************************************
// Try to implement this in fortran using C bindings       *
// or use standard fortran routines directly if available. *
// We need only these routines:                            *
// system: why is it needed? why call system is not OK     *
// getpwnam()->pw_dir - home directory from /etc/passwd    *
// getpid()                                                *
// getenv()                                                *
// putenv()                                                *
// getpwuid(),getlogin() - get username                    *
// OUT::etime() - use standard fortran                     *
// time(),ctime() for fdate - use fortran!                 *
// system() for csystem ??                                 *
// Not needed::: gettimeofday() for xsecnds - used!        *
// uname(),gethostbyname() used for call uninf! we need!   *
// - what to do with test endian (can we do it different?) *
// - read namd: maybe replace by access='stream'           *
//                                                         *
// *********************************************************

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// begin for get_username
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
// end for get_username

// for get_home_dir: wordexp, tolower
#include <wordexp.h>
#include <ctype.h>

// for backtrace
#include <execinfo.h>

// aag 06/07
// read/write binary namd file
void readnamd(double * x, double * y, double * z, int * ptr_natom,
              char * fname, int * ptr_flen, int * ptr_ier) {
  FILE * fp;
  char * filename;

  // should not be INT since this can actually be long
  //   not OK in binary file
  int n, i;

  // allocate space for file names
  if ( (filename = malloc(*ptr_flen + 1)) == NULL) {
    printf("Error from malloc in readnamd\n");
    *ptr_ier = -1;
    return;
  }

  // now copy the strings and null terminate them
  strncpy(filename,fname,*ptr_flen);
  filename[*ptr_flen] = '\0';
  printf("opening file %s\n",filename); 
  if ( (fp = fopen(filename,"rb")) == NULL) {
    printf("could not open file %s\n",filename); 
    free(filename);
    *ptr_ier = -1;
    return;      
  }

  fread(&n, sizeof(int), 1, fp);

  // checking if the number of atoms is the same
  if (n != (*ptr_natom)) {
    printf("number of atoms does not match; in file %d, in psf %d ; exiting...\n",
          (int) n, (int) *ptr_natom);
    free(filename);
    *ptr_ier = -1;
    return;      
  }

  // read in coordinates
  for (i = 0; i < n; i++) {
    fread(x + i, sizeof(double), 1, fp);
    fread(y + i, sizeof(double), 1, fp);
    fread(z + i, sizeof(double), 1, fp);
  }
  fclose(fp);

  free(filename);
  *ptr_ier = 0;

}

void writenamd(double * x, double * y, double * z, int * ptr_natom,
               char * fname, int * ptr_flen, int * ptr_ier)
{
  FILE * fp;
  int n, i;
  char * filename;
  
  // allocate space for file names
  if ( (filename = malloc(*ptr_flen + 1)) == NULL) {
    printf("Error from malloc in writenamd\n");
    *ptr_ier = -1;
    return;
  }
  
  // now copy the strings and null terminate them
  strncpy(filename,fname,*ptr_flen);
  filename[*ptr_flen] = '\0';
  printf("opening file %s\n",filename); 
  if ( (fp = fopen(filename,"wb")) == NULL) {
    printf("could not open file %s\n",filename); 
    free(filename);
    *ptr_ier = -1;
    return;      
  }

  n = *ptr_natom;
  printf("saving number of atoms = %d\n", (int) n);
  fwrite(&n, sizeof(int), 1, fp);
  
  // write in coordinates
  for (i = 0; i < n; i++) {
    fwrite(x + i, sizeof(double), 1, fp);
    fwrite(y + i, sizeof(double), 1, fp);
    fwrite(z + i, sizeof(double), 1, fp);
  }
  fclose(fp);

  free(filename);
  *ptr_ier = 0;

}

int unbuffer_stdout() { 
  return setvbuf(stdout, NULL, _IONBF, 0);
}

int get_username(char * name, int * namelen) {
  if (*namelen <= 0) return -1;
  name[0] = '\0';
  
  char * env_name = getenv("USER");
  if (env_name != NULL && env_name[0] != '\0') {
    strncpy(name, env_name, *namelen);
    *namelen = strnlen(name, *namelen);
  }
#if STATIC == 1
  else {
    *namelen = 0;
    return -1;
  }
#else
  else {
    struct passwd *pws;
    pws = getpwuid(geteuid());

    if (pws != NULL && pws->pw_name[0] != '\0') {
      strncpy(name, pws->pw_name, *namelen);
      *namelen = strnlen(name, *namelen);
    } else {
      *namelen = 0;
      return -1;
    }
  }
#endif
  return 0;
}

int expand_tilde(char * in_exp, int in_length,
                 char * out_exp, int * out_length) {
  // Initialize result in case of error
  if (*out_length > 0) {
    out_exp[0] = '\0';
  }

  if (*out_length <= 0 || in_length <= 0 ||
      in_exp == NULL || in_exp[0] == '\0'
      || out_exp == NULL) {
    *out_length = 0;
    return -1;
  }

  // make sure wordexp input is properly null terminated
  // because this function will be called from fortran
  int safe_length = strnlen(in_exp, in_length);
  if (safe_length <= 0) {
    *out_length = 0;
    return -1;
  }

  char * safe_exp = (char *) malloc((1 + safe_length) * sizeof(char));
  strncpy(safe_exp, in_exp, safe_length);
  
  safe_exp[safe_length] = '\0';
  safe_length = strnlen(safe_exp, safe_length);

  if (safe_length < in_length) {
    *out_length = 0;
    return -1;
  }

  wordexp_t exp_result;
  int status = wordexp(safe_exp, &exp_result, 0);
  free(safe_exp);

  if (status != 0 || exp_result.we_wordc <= 0 || exp_result.we_wordv[0] == NULL) {
    *out_length = 0;
    wordfree(&exp_result);
    return -1;
  }

  strncpy(out_exp, exp_result.we_wordv[0], *out_length);
  wordfree(&exp_result);

  *out_length = strnlen(out_exp, *out_length);
  return 0;
}

// This function is better than system() because it does not
// replicate program's memory requirements
int fsystem(char * com, int com_length)
{
  FILE * a;
  char c;
  int stat;

  com[com_length] = 0x0;
  if ( (a = popen(com, "r")) == NULL ) return -1;
  
  // Read the shell's stream and put it to the stdout
  while ( (c = getc(a) ) != EOF ) putchar(c) ;

  // Wait to finish with this command
  stat = pclose(a) ;

  // flush the stdout buffer to see the results
  // immediately instead of waiting till the end
  // of the CHARMM script
  fflush(stdout);
  
  return(stat);
}

void fputenv(char * env, int len) {
  int err;
  size_t p;
    
  char * pp = strchr(env, '=');
  if (pp == NULL) {
    perror("fputenv cannot find = in env statement");
    exit(110);
  }
  
  p = pp - env;

  char * name = (char *) malloc(len * sizeof(char));
  if (name == NULL) {
    perror("fputenv cannot allocate memory");
    exit(110);
  }
  
  strncat(name, env, p);
  name[p] = '\0';
  
  char * val = (char *) malloc(len * sizeof(char));  
  if (val == NULL) {
    perror("fputenv cannot allocate memory");
    exit(110);
  }

  strncpy(val, &env[p + 1], len);
  
  err = setenv(name, val, 1);

  free(name);
  free(val);
  
  if (err != 0 ) {
    printf("fputenv=%d\n", err);
    perror("fputenv problem");
    exit(110);
  }
}

void stack_trace() {
  void * callstack[128];
  int i, frames = backtrace(callstack, 128);
  char ** strs = backtrace_symbols(callstack, frames);

  for (i = 0; i < frames; ++i) {
    printf("%s\n", strs[i]);
  }
  free(strs);
}
