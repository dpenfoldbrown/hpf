/** Show cluster information for decoys in a silent mode file; **/

#include <stdio.h>
#include <stdlib.h>
#include "rms2.c"
#include <assert.h>
#include <string.h>

#define MAX_NUM_DECOYS  100000
#define MAX_NUM_CLUSTERS   1000
#define MAX_COMMANDFILE_LINE_LENGTH 10000

/*  #define DEFAULT_MIN_THRESHOLD 4.0  NOW WE ARE USING RMS100default_min*/
#define DEFAULT_MIN_THRESHOLD_RMS100 5.0
#define DEFAULT_MAX_THRESHOLD 12.0

#define MAX_SEQUENCE_LENGTH 1000
#define MAX_NUM_HOMOLOGS    200


int Mcount = 0;

/*******************************************************/
/** GLOBAL VARIABLES (probably bad form, but oh well) **/

double **commonCoords, **fullCoords, **decoyScores, *rmsdToNative,
  *goodDecoyCoords, /* for doing rmsds to native */
  *goodNativeCoords; /* stores native coords for (possibly a subset of) common positions */

char **commonStructure, **decoyName, **homologName, **homologSequence, *ssConstraint;
int *homologLength, **commonToFull, **fullToCommon, *decoyToHomolog,
  *goodCoords, /* remembers which positions in commonCoords have native coordinates */
  *commonToGood, /* remembers which index in the good coords set: like commonToGood */
  argCount;

char *outputFilename; /** filename for output, get it from commandfile **/

int numResidues,numHomologs,NATIVE,numDecoys,numGoodCoords,
  minClusterSize,minTopClusterSize,targetClusterSize,maxClusterSize;



float minThreshold,maxThreshold;


/** structures **/
struct IntList {
  int              data;
  struct IntList  *next;
};

typedef struct IntList    INT_LIST;
typedef INT_LIST         *INT_LINK;



/** functions **/

int IntMin(int a, int b) { 
  return (a<b?a:b);
}

float RMS100(float rmsd,int L) {
  if (L<14) return 1000;
  else return rmsd/ (1 + 0.5 * (log (( (float) L)/ 100.0)));
}

int SS2Int(char a) {

  if (a=='E') return 0;
  else if (a=='H') return 1;
  else if (a=='L') return 2;
  else {
    fprintf(stderr,"?2L");
    return 2;
  }
}

int CompareIncreasing(float **a,float **b) {
  if (a[0][0]>=b[0][0]) return 1;
  else return -1;
}
int CompareIntIncreasing(int *a,int *b) {
  if (a[0]>=b[0]) return 1;
  else return -1;
}

int CompareDecreasing(float **a,float **b) {
  if (a[0][0]>=b[0][0]) return -1;
  else return 1;
}

int IntListLength(INT_LINK l) {
  if (l==NULL) return 0;
  return (1+IntListLength(l->next));
}

void DeleteIntList(INT_LINK l) {
  if (l==NULL) return;
  DeleteIntList(l->next);
  free(l);
  return;
}

INT_LINK AddIntLink(int i,INT_LINK l) {
  INT_LINK new;

  new = malloc(sizeof(INT_LIST));
  Mcount ++;
  new->next = l;
  new->data = i;
  return new;
}




void ReadSilentFile(int homolog) {
  FILE *file;
  char junk50[10000],*ss,*filename,sequence[MAX_SEQUENCE_LENGTH];
  int i,j,L,rsdCount,startNumDecoys,totalDecoys;

  
  filename = homologName[homolog];
  startNumDecoys = numDecoys;

  if ( (file = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"cant open file: %s\n",filename);
    return;
  }
  
  fscanf(file,"%s %s%*c",junk50,sequence);
  L = strlen(sequence);

  assert (!strcmp (homologSequence[homolog], sequence));

  fprintf(stderr,"LOAD  %s L = %d ",
	  homologName[homolog],L);
  
  /** read scoreline header **/

  fseek(file,0,0);
  fscanf(file,"%*[^\n]%*c%*[^\n]%*c");
  
  ss = calloc(L+5,sizeof(char));

  while (fscanf(file,"%s",junk50)==1) {
    if (strlen(junk50) > 5000) fprintf(stderr,"whoah: %s\n",junk50);
    if (strcmp(junk50,"SCORE:")!=0) continue;

    fullCoords[numDecoys] = calloc(L*3,sizeof(double));
    
    fscanf(file,"%*[^\n]%*c");
    
    rsdCount = 0;
    while (fscanf(file,"%d %c",&i,ss+rsdCount)==2) {
      if ( !(ss[rsdCount] == 'H' || ss[rsdCount] == 'L' || ss[rsdCount] == 'E')) break;
      if (i!=rsdCount+1 || i> L) break;
      if (fscanf(file,"%*f %*f %*f %lf %lf %lf %[^\n]%*c",
		 fullCoords[numDecoys]+3*rsdCount,
		 fullCoords[numDecoys]+3*rsdCount+1,
		 fullCoords[numDecoys]+3*rsdCount+2,junk50)!=4) break;
      if (strlen(junk50) > 5000) fprintf(stderr,"whoah: %s\n",junk50);
      rsdCount++;
    }
    
    if (i!=L || rsdCount!=L) {
      /** Some problem **/
      fprintf(stderr,"*(%s)",junk50);
/*        fprintf(stderr,"silentfile error: %d %d %d %s\n",i,rsdCount,L,junk50);  */
      free(fullCoords[numDecoys]);
      continue;
    }

    ss[L] = '\0';
      
    /** check secondary structure constraint **/
    assert (strlen(ss)==L);

    j = strlen(junk50);
    for (i=0;i<j;i++) /** Convert name to alpha-numeric+_ (fill in spaces,etc) **/
      if ( ! ( (junk50[i]>='A' && junk50[i]<='Z') ||
	       (junk50[i]>='a' && junk50[i]<='z') ||
	       (junk50[i]>='0' && junk50[i]<='9') || 
	       (junk50[i] == '.' ||
		junk50[i] == '/'))) 
	junk50[i] = '_';
    
    commonStructure[numDecoys] = calloc(numResidues+1,sizeof(char));

    /** assign sec structure **/
    for (i=0;i<numResidues;i++) 
      commonStructure[numDecoys][i] = ss[ commonToFull [homolog] [i] ];
    commonStructure[numDecoys][numResidues] = '\0';

    /** assign decoy name **/
    decoyName[numDecoys] = calloc( 10 + strlen(junk50),sizeof(char));
    sprintf(decoyName[numDecoys],"%d,%s",numDecoys+1,junk50);

/*      decoyName[numDecoys] = calloc( 10 + strlen(junk50) + strlen(homologName[homolog])+5, */
/*  				   sizeof(char)); */
/*      sprintf(decoyName[numDecoys],"%d,%s,%s",numDecoys+1,homologName[homolog],junk50); */

    /** warn about possible wacked coordinates **/ 
    if (fullCoords[numDecoys] [3*L-3] == 0.0 && 
	fullCoords[numDecoys] [3*L-2] == 0.0 &&
	fullCoords[numDecoys] [3*L-1] == 0.0 ) {
      fprintf(stderr,"\nWARNING: %s:%s coordinates end in 0.000's -- possible error?\n",
	      homologName[homolog],junk50);
    }

    /** assign commonCoords **/
    commonCoords[numDecoys] = calloc(numResidues*3,sizeof(double));
    for (i=0;i<numResidues;i++) 
      for (j=0;j<3;j++)
	commonCoords[numDecoys][3*i+j] = fullCoords[numDecoys][3 * (commonToFull[homolog][i])+j];
  
    decoyToHomolog[numDecoys] = homolog;

    /** progress report **/
    if ( (numDecoys%100)==0) fprintf(stderr,"-");

    numDecoys++;
    if (numDecoys >= MAX_NUM_DECOYS) {
      fprintf(stderr,"STOP: numDecoys >= MAX_NUM_DECOYS\n");
      break;
    }
  } /** while(fscanf... **/
  
  
  totalDecoys = numDecoys - startNumDecoys; /** # decoys successfully read from this file **/
  fprintf(stderr," total: %d %s\n",
	  totalDecoys,homologName[homolog]);
  
  fclose(file);
  
  return;
}

void FreeNeighbors (float ***neighbors,int N, int L) {
  int i,j;

  for (i=0;i<N;i++) {
    for (j=0;j<L;j++) 
      free(neighbors[i][j]);
    free(neighbors[i]);
  }
}

float *** MakeNeighborList (int N,double **coords, int L, int start, int stop,
			    int monitorProgress) { 
  /* N=#coords, L=#neighbors */ 

  float ***neighbors,**nl,r;
  int f1,f2,I,J,i,j,numPositions;
  
  numPositions = stop-start+1;

/*    fprintf(stderr,"Computing decoy-decoy rmsds: numDecoys=%d, numNeighbors=%d, numPositions=%d.\n", */
/*  	  N,L,numPositions); */
  
  /* initialize the list: */
  neighbors = calloc(N,sizeof(float**));
  for (f1=0;f1<N;f1++) {
    /** neighbors[f1] is a sorted list of [rmsd,f2] to the nearest L nbrs seen **/
    neighbors[f1] = calloc(L,sizeof(float*));
    for (j=0;j<L;j++) {
      neighbors[f1][j] = calloc(2,sizeof(float));
      neighbors[f1][j][0] = 1000.0;
      neighbors[f1][j][1] = -1;
    }
  }
  if (monitorProgress) {
    for (f1=0;f1<N;f1++) 
      if (  !(f1%10) ) fprintf(stderr,"-");
    fprintf(stderr,"\n");
  }
  
  /* calculate rmsd's */
  for (f1=0;f1<N;f1++) {
    
    if ( monitorProgress && !(f1%10) ) fprintf(stderr,".");

    for (f2 = f1+1;f2<N;f2++) {
	
      r = (float) rmsfit_(&numPositions, coords[f1]+3*start, coords[f2]+3*start);
      for (I=f1,J=f2;
	   I<=f2;
	   I+=(f2-f1),J-=(f2-f1)) {

	nl = neighbors[I];
	
	if (r < nl[L-1][0]) {
	  nl[L-1][0] = r;
	  nl[L-1][1] = J;
	  i = L-1;
	  while (i>0 && nl[i-1][0] > r) {
	    nl[i][0] = nl[i-1][0];
	    nl[i][1] = nl[i-1][1];
	    nl[i-1][0] = r;
	    nl[i-1][1] = J;
	    i--;
	  }
	}
      }
    }
  }

  if (monitorProgress) fprintf(stderr,"\n");
  return neighbors;
}


float ClusteringThreshold(float ***neighbors, int N,
			  int minCS, int targetCS, int maxCS,
			  float minT, float maxT,
			  int *clusterCenterPointer, 
			  int *targetPointer,
			  int monitorProgress) {

  int f1,clusterCenter, target, lastTarget;
  float r, bestRMSD, threshold;

  target = targetCS;
  lastTarget = 0;

  if (monitorProgress) fprintf(stderr,"Find threshold: (target,threshold) ");
  while (target <= maxCS && target >= minCS) {
    bestRMSD = 100.0;
    for (f1=0;f1<N;f1++) {
      r = neighbors[f1][target-1][0]; /** the rmsd of (target+1)th neighbor. **/
      if (r<bestRMSD) {
	bestRMSD = r;
	clusterCenter = f1;
      }
    }
    if (monitorProgress) fprintf(stderr," %d,%f",target,bestRMSD);
    
    if (bestRMSD < minT) {
      if (lastTarget == target+1) break; /** bouncing backward **/
      lastTarget = target;
      target++;
    }
    else if (bestRMSD > maxT) {
      if (lastTarget == target-1) break; /** bouncing backward **/
      lastTarget = target;
      target--;
    }
    else break;
  }
  if (monitorProgress) fprintf(stderr,"\n");
  
  if (target>maxCS) target = maxCS;
  else if (target<minCS) target = minCS;
  
  threshold = bestRMSD;
  assert (target<= maxCS && target >= minCS);

  *targetPointer = target;
  *clusterCenterPointer = clusterCenter;
  
  return threshold;
}

/** returns 1 if we cant open the file **/
integer ReadCommandFile(char *filename) {
  
  int i,j,homologCount,foundTarget,alignmentLength,homolog,length,commonIndex,
    fullIndex,common[MAX_COMMANDFILE_LINE_LENGTH],MIN_THRESHOLD_SET_IN_COMMAND_FILE;
  
  float minThresholdRMS100;
  
  char line[MAX_COMMANDFILE_LINE_LENGTH], *alignment[MAX_NUM_HOMOLOGS],
    token[MAX_COMMANDFILE_LINE_LENGTH];
  
  FILE *file;

  if ( (file=fopen(filename,"r"))==NULL) {
    fprintf(stderr,"cant open command file: %s\n",filename);
    return 1;
  }
  
  fprintf(stderr,"Reading command file: %s\n",filename);

  for (i=0;i<MAX_COMMANDFILE_LINE_LENGTH;i++) common[i] = 1;
  
  homologCount = 0;
  foundTarget = 0;
  alignmentLength = 0;

  /** user can set these parameters **/
  maxThreshold = DEFAULT_MAX_THRESHOLD;
  MIN_THRESHOLD_SET_IN_COMMAND_FILE = 0; /** if ==0 after read file use rms100 **/
  minThresholdRMS100 = DEFAULT_MIN_THRESHOLD_RMS100;
  

  outputFilename = calloc(50,sizeof(char));
  sprintf(outputFilename,"./cluster_output.txt");
  
  
  while (fscanf(file,"%[^\n]%*c",line)==1) { /** read a single line **/
    if (strlen(line) > MAX_COMMANDFILE_LINE_LENGTH-5) {
      fprintf(stderr,"commandfile line too long; recompile robetta_cluster\n");
      fclose(file);
      return 1;
    }
    
    sscanf(line,"%s",token);
    
    if (!strcmp(token,"MIN_THRESHOLD")) {
      MIN_THRESHOLD_SET_IN_COMMAND_FILE=1;
      sscanf(line,"%*s %f",&minThreshold);
    }

    if (!strcmp(token,"MIN_THRESHOLD_RMS100")) {
      sscanf(line,"%*s %f",&minThresholdRMS100);
    }
    
    if (!strcmp(token,"MAX_THRESHOLD")) {
      sscanf(line,"%*s %f",&maxThreshold);
    }
    
    if (!strcmp(token,"OUTPUT_FILE")) {
      sscanf(line,"%*s %s",token);
      outputFilename = calloc(strlen(token)+5,sizeof(char));
      strcpy(outputFilename,token);
    }
    
    if ((!strcmp(token,"TARGET")) || (!strcmp(token,"HOMOLOG"))) {
      
      
      if (!strcmp(token,"TARGET")) {
	foundTarget = 1;
	homolog = 0;
      }
      else {
	homologCount = homologCount+1;
	homolog = homologCount;
      }

      sscanf(line,"%*s %s",token); /** the filename **/

      
      homologName[homolog] = calloc(strlen(token)+5,sizeof(char));
      strcpy(homologName[homolog],token);

      sscanf(line,"%*s %*s %s",token); /** the alignment **/
      
      if (alignmentLength) assert (strlen(token) == alignmentLength);
      else alignmentLength = strlen(token);
      
      alignment[homolog] = calloc(alignmentLength+1,sizeof(char));
      strcpy(alignment[homolog],token);
      
      length = 0;
      homologSequence[homolog] = calloc(alignmentLength+1,sizeof(char));
      for (j=0;j<alignmentLength;j++) {
	if (token[j] == '-' || token[j] == '.') 
	  common[j] = 0;
	else {
	  homologSequence[homolog][length] = token[j];
	  length++;
	}
      }
      homologSequence[homolog][length] = '\0';
      homologLength[homolog] = length;
    }
  }
  
  if (!foundTarget) {
    fprintf(stderr,"Couldnt find TARGET line in commandfile: %s\n",filename);
    return 1;
  }
  
  numResidues = 0;
  for (i=0;i<alignmentLength;i++) numResidues+=common[i];

  numHomologs = homologCount+1;
  
  for (i=0; i<numHomologs; i++) {
    commonToFull[i] = calloc(alignmentLength,sizeof(int));
    fullToCommon[i] = calloc(alignmentLength,sizeof(int));

    commonIndex = 0;
    fullIndex = 0;
    for (j=0;j<alignmentLength;j++) {
/*        fprintf(stderr,"%d %d\n",i,j); */
      if (alignment[i][j] == '-' || alignment[i][j] == '.') continue;
      else if (common[j]) {
	commonToFull[i][commonIndex] = fullIndex;
	fullToCommon[i][fullIndex] = commonIndex;
/*  	fprintf(stderr,"%d %d %d-%d\n",i,j,fullIndex,commonIndex); */
	commonIndex++;
	fullIndex++;
      }
      else {
	fullToCommon[i][fullIndex] = -1;
/*  	fprintf(stderr,"%d %d %d-\n",i,j,fullIndex); */
	fullIndex++;
      }
    }
    assert (fullIndex == homologLength[i] && commonIndex == numResidues);
  }

  for (i=0;i<numHomologs;i++) {
    fprintf(stderr,"%3d %s %s\n",i,alignment[i],homologName[i]);
  }
  fprintf(stderr,"cmn:");
  for (i=0;i<alignmentLength;i++) {
    if (common[i]) fprintf(stderr,"*");
    else fprintf(stderr,"-");
  }
  fprintf(stderr,"\n");
    
  numDecoys = 0;
  
  /** Read the silent files **/
  for (i=0;i<numHomologs;i++) {

    ReadSilentFile(i);
  }


  /** setup the parameters **/
  
  /* minThreshold is set in command file; default is by rms100, see below **/
  /* maxThreshold is set in command file; default is DEFAULT_MAX_THRESHOLD **/
  
  minClusterSize = 5;
  minTopClusterSize = 20;
  targetClusterSize = numDecoys / 50; /** 40  for 2000,  80 for 4000 **/
  maxClusterSize = (3*numDecoys)/40; /**  150 for 2000, 300 for 4000 **/
  
  /** do some checks, in case numDecoys is small **/
  if (targetClusterSize < minTopClusterSize) targetClusterSize = minTopClusterSize;
  if (maxClusterSize < 2*minTopClusterSize) maxClusterSize = 40;
  
  
  if (! MIN_THRESHOLD_SET_IN_COMMAND_FILE) { /** use rms100 to set minThreshold **/
    
    minThreshold = minThresholdRMS100 * 
      (1.0 + 0.5 * log(( (float) numResidues)/100.0));
    
    if (minThreshold < 3) minThreshold = 3;

    fprintf(stderr,"Set minThreshold = %f = rms100 of %f, based on %d residues\n",
	    minThreshold,
	    minThresholdRMS100,
	    numResidues);
    
  }
  
  if (maxThreshold < minThreshold + 2) {
    maxThreshold = minThreshold + 2;
    fprintf(stderr,"WARNING: maxThreshold is being reset to minThreshold + 2\n");
  }

  NATIVE = 0; /** currently, not configured to use native info **/

  fprintf(stderr,"robetta_cluster parameter: minThreshold= %f\n",minThreshold);
  fprintf(stderr,"robetta_cluster parameter: maxThreshold= %f\n",maxThreshold);
  fprintf(stderr,"robetta_cluster parameter: minTopClusterSize= %d\n",minTopClusterSize);
  fprintf(stderr,"robetta_cluster parameter: minClusterSize= %d\n",minClusterSize);
  fprintf(stderr,"robetta_cluster parameter: targetClusterSize= %d\n",targetClusterSize);
  fprintf(stderr,"robetta_cluster parameter: maxClusterSize= %d\n",maxClusterSize);
  fprintf(stderr,"robetta_cluster parameter: outputFilename= %s\n",outputFilename);
  
  
  
  return 0;
}




/** returns the lists of cluster members: first member in each list is size, 2nd is center **/

int ** DoClustering( float ***neighbors, float threshold, int target, int clusterCenter, 
		     int* numClustersPointer) {

  int i,j,I,f1,clusterCount,count;

  int *exists, **clusterMembers;
  float **clusterList,**nl;

  /** setup **/
  exists = calloc(numDecoys,sizeof(int)); /** remembers who is already in a cluster **/
  for (i=0;i<numDecoys;i++) exists[i] = 1;
  
  clusterList = calloc(numDecoys,sizeof(float*)); /** for sorting centers by cluster size **/
  for (i=0;i<numDecoys;i++) {
    clusterList[i] = calloc(2,sizeof(float));
  }

  clusterMembers = calloc(MAX_NUM_CLUSTERS,sizeof(int*)); /** list of all cluster members **/
  for (i=0;i<MAX_NUM_CLUSTERS;i++)
    clusterMembers[i] = calloc(target+5,sizeof(int));
  
  /** testing **/
  for (i=0;i<numDecoys;i++) {
    for (j=0;j<maxClusterSize;j++) {
      I = (int) neighbors[i][j][1];
      assert ( fabs(neighbors[i][j][1]-I) < 0.001 );
    }
  }
  
  clusterCount = 0;


  while (1) { /** make clusters in this loop -- break when size is < minClusterSize **/
    
    if (clusterCount) {
      clusterCenter = (int) clusterList[0][1];
    }

    
    clusterMembers[clusterCount][1] = clusterCenter;
    exists[clusterCenter] = 0;

    count = 1;
    nl = neighbors[clusterCenter];
    for (i=0;i<target;i++) {
      if (nl[i][0] >= threshold) break;
      if (! exists[ (int) nl[i][1]]) continue;
      clusterMembers[clusterCount][count+1] =  (int) nl[i][1];
      exists[ (int) nl[i][1] ] = 0;
      count++;
    }
    
    clusterMembers[clusterCount][0] = count; /** first element in list is size **/
    
    /** sanity check **/
    if (clusterCount) {
      assert (count==clusterList[0][0]);
    }
    else if (count!=target) fprintf(stderr,"WARNING: count != target\n");

    
    /** Update clusterList **/
    for (f1=0;f1<numDecoys;f1++) {
      

      clusterList[f1][1] = f1;
      if (! exists[f1] ) {
	clusterList[f1][0] = 0;
	continue;
      }

      count = 1;
      nl = neighbors[f1];
      for (i=0;i<target;i++) {
	if (nl[i][0] >= threshold) break;
	if (exists[ (int) nl[i][1]])
	  count++;
      }

      clusterList[f1][0] = count;
    }
    
    qsort(clusterList,numDecoys,sizeof(float*),
	  (int (*)(const void*,const void*))CompareDecreasing);
    
    clusterCount++;

    if (clusterList[0][0] < minClusterSize) break;
    
    if (clusterCount>=MAX_NUM_CLUSTERS) {
      fprintf(stderr,"Clustering was stopped at %d clusters\n",MAX_NUM_CLUSTERS);
      break;
    }
  }

  *numClustersPointer = clusterCount;
  return clusterMembers;
}
  
void MakeOutput (char *filename, float threshold, int target, 
		 int **clusterMembers, int numClusters) {


  int i,cluster,*members,size;
  FILE *file;

  if ( (file=fopen(filename,"w"))==NULL) {
    fprintf(stderr,"Couldnt open output file: %s\n",filename);
    return;
  }
  
  fprintf(file,"THRESHOLD: %f TOP_CLUSTER_SIZE: %d\n",threshold,target);

  fprintf(file,"#: center size cluster-cluster-rmsds\n");
  
  for (cluster=0;cluster<numClusters;cluster++) {

    members = clusterMembers[cluster];
    size = members[0];
    members++;
    
    fprintf(file,"%d: %s %d",cluster,decoyName[members[0]],size);
    for (i=0;i<=cluster;i++) 
      fprintf(file," %9.3f",rmsfit_(&numResidues,
					commonCoords[members[0]],
					commonCoords[clusterMembers[i][1]]));
    fprintf(file,"\n");
  }

  
  fprintf(file,"\n\nCLUSTER MEMBERS:\n");
  for (cluster=0;cluster<numClusters;cluster++) {

    members = clusterMembers[cluster];
    size = members[0];
    members++;
    fprintf(file,"%d:",cluster);
    for (i=0;i<size;i++) 
      fprintf(file," %s",decoyName[members[i]]); /** numbering starts at 0**/
    fprintf(file,"\n");
  }
  fclose(file);
}
      
	    

/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/


int main(int argc, char *argv[]) {

  float ***neighbors,threshold;

  int **clusterMembers,done;

  int clusterCenter,numClusters,
    target;


  
  
  /** a bit of setup **/
  homologName = calloc(MAX_NUM_HOMOLOGS,sizeof(char*));
  homologLength = calloc(MAX_NUM_HOMOLOGS,sizeof(int));
  homologSequence = calloc(MAX_NUM_HOMOLOGS,sizeof(char*));
  commonToFull = calloc(MAX_NUM_HOMOLOGS,sizeof(int*));
  fullToCommon = calloc(MAX_NUM_HOMOLOGS,sizeof(int*));
  homologName++; /** native index is -1 **/
  homologLength++;
  homologSequence++;
  commonToFull++;
  fullToCommon++;



  commonCoords = calloc(MAX_NUM_DECOYS,sizeof(double*));
  commonStructure = calloc(MAX_NUM_DECOYS,sizeof(char*));
  fullCoords = calloc(MAX_NUM_DECOYS,sizeof(double*));
  decoyScores = calloc(MAX_NUM_DECOYS,sizeof(double*));
  decoyName = calloc(MAX_NUM_DECOYS,sizeof(char*));
  decoyToHomolog = calloc(MAX_NUM_DECOYS,sizeof(int));

  commonCoords++;  /** native index is -1 **/
  commonStructure++;
  fullCoords++;
  
  decoyToHomolog[-1] = -1;



  if (ReadCommandFile( argv[1] )) return 1;
  


  /****************************************************************************/
  /*********************************** clustering *****************************/
  
  
  /** calculate decoy-decoy rmsds  ---- THIS IS THE SLOW STEP **/
  fprintf(stderr,"decoy-decoy RMSDS\n");


  done = 0;

  while (!done) { /** keep looping until we have threshold >= minThreshold, probably only once*/
    neighbors = MakeNeighborList (numDecoys, commonCoords, maxClusterSize,0,numResidues-1,1);
  

    /** find clustering threshold **/
    threshold = ClusteringThreshold (neighbors, numDecoys, 
				     minTopClusterSize, targetClusterSize, maxClusterSize,
				     minThreshold, maxThreshold,
				     &clusterCenter, 
				     &target,
				     1); /** show progress reports **/
    if ( (threshold < minThreshold) && (maxClusterSize<numDecoys-1) ) {
      FreeNeighbors(neighbors,numDecoys,maxClusterSize);

      /** we set maxClusterSize too small **/
      maxClusterSize += (numDecoys/10);
      maxClusterSize = (maxClusterSize > numDecoys-1)?(numDecoys-1):maxClusterSize;
      
      fprintf(stderr,
	      "threshold= %f minThreshold= %f repeat clustering with new maxClusterSize= %d\n",
	      threshold, minThreshold,maxClusterSize);
    }
    else done = 1;
  }
  
  clusterMembers = DoClustering (neighbors, threshold, target, clusterCenter,&numClusters);
  
  assert (clusterMembers[0][1] == clusterCenter); /* first member is center of each cluster*/
  
  
  MakeOutput (outputFilename, threshold, target, clusterMembers, numClusters);



  fprintf(stderr,"clustering_threshold: %f top_cluster_size: %d top_cluster_center: %d \n",
	  threshold,target,clusterCenter);
  
  return 0;
}








