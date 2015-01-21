#include <sys/types.h>
#include <sys/wait.h>
#include <signal.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <fcntl.h>


int main(int argc, char* argv[]) {


 
  if (argc < 2) {
    printf("usage: %s <command> <timout(s)> \n", argv[0]);
    exit(EXIT_FAILURE);
  }

  //Locating the binary file for dymosim
  char *path;
  path = getenv ( "PWD" );
  char dymcom[1024];
  snprintf(dymcom,sizeof(dymcom),"%s/dymosim",path); 
  
  
  //Initializing the arguments for the dymola child 
  

  char **dym; // = {"dymosim ",dymsimarg,(char *)0};
  dym = (char **) malloc(sizeof(char *) * argc); 
  dym[0] = "dymosim"; 
  dym[argc-1] = (char *) 0;
  int i; 
  for(i=1;i<argc-1;i++)
  dym[i] = argv[i]; 

 
  unsigned long sleeptime    =   strtoul(argv[argc-1], NULL, 10);

  if (sleeptime == 0) {
    printf("Invalid timeout argument\n");
    printf("usage: %s <dymosimarg> <timout(s)> \n", argv[0]);
    exit(EXIT_FAILURE);
  }


  //This is not a child of any body then !? 
  /*int outfd = open("dymosim.out", O_CREAT | O_WRONLY, 00664);
  if (outfd == -1) {
     printf("Failed to open dymosim.out\n");
     exit(EXIT_FAILURE);
  }
  dup2(outfd, STDOUT_FILENO);
  dup2(outfd, STDERR_FILENO);*/

 
  pid_t p = fork();
  pid_t p2 = 10;  



  if (p != 0) {
    
    printf("[Parant] Dymola Child id: %u\n",p);
    p2 = fork(); 
    
    if(p2 != 0) { 
	 
      printf("[Parant] Sleepy Child id: %u\n",p2);
      int status = 0;
      // Wait any child process 
#ifdef DEBUG
       printf("[parent] waitpid()\n");
#endif
      pid_t rv = waitpid(-1,&status,0);

#ifdef DEBUG
      printf("[Parant] Waiting : %d\n",rv);
#endif

      if(rv < 0) { 

      } else { // See which process has finished and kill the other 
	int kstatus; 

	if(rv == p) { 
	  kstatus = kill(p2,SIGKILL);//SIGHUP);
#ifdef DEBUG
	  printf("[Parant] Sleepy Child got killed :%d\n",kstatus); 
#endif
	} else if (rv == p2) {
	  kstatus = kill(p,SIGKILL);
	    
#ifdef DEBUG
	  printf("[Parant] Dymola Child got killed :%d\n",kstatus);
#endif   
	}
      } // if Kill

    } else { // Create a sleepy child 

#ifdef DEBUG
      printf("[sleepy Child] Sleeping %u second\n", (unsigned) sleeptime);
#endif
      sleep(sleeptime);
	

#ifdef DEBUG
      printf("[Sleepy Child]: I am finished\n");
#endif
      exit(EXIT_SUCCESS);
	 
    } // if(p2 != 0) 
      
   
  }  // if p == 0 Dymola child 
  else {

#ifdef DEBUG
    printf("[Dymola Child] : Starting work\n");
#endif
       
    int  execstatus;
    execstatus = execv(dymcom,dym);
    
#ifdef DEBUG
    printf("[Dymola Child] Dymola command : %d\n",execstatus);
#endif
   
    exit(EXIT_SUCCESS);
  }

  //close(outfd);

  free(dym);
 
#ifdef DEBUG
  printf("Thread: %u %u\n",p,p2);
#endif 
  exit(EXIT_SUCCESS);

}


