// C program to input an ascii picture of a tangle
// and output the matrix it would compute in a Fibonacci quantum computer.
// So far it can handle crossings and cups and caps.
//
// We work over the field of integers mod 521, with q=5,
// chosen because we need q^5=-1 and q+q^-1 has a square root.
// Basis vectors are enumerated by binary numbers with no consecutive 1's.
// 
// Input:
// |, \, / mean a vertical strand
// % means a positive crossing
// S means a negative crossing
// u means a cup
// n means a cap
// , is the same as "new line"
// . means end (the input should end with ",.")
// spaces are ignored
 
// Example input:  |% , %| , |% , .

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MOD 521
// Modulo 521
// q=5
// [2] = 422
// sqrt[2] = 237

// The fibonacci matrix for a crossing is:
// 1/[2]        0    q^8/sqrt[2]  0    0
// 0            q^2  0            0    0
// q^8/sqrt[2]  0    q/[2]        0    0
// 0            0    0            q^2  0
// 0            0    0            0    q^9
// 
#define CROSS {421,0,94,0,0, 0,25,0,0,0, 94,0,21,0,0, 0,0,0,25,0, 0,0,0,0,417}
#define UNCROSS {421,0,398,0,0, 0,396,0,0,0, 398,0,501,0,0, 0,0,0,396,0, 0,0,0,0,5}

// The matrix for a cup is:
// 1         0
// 0         0
// 1/sqrt[2] 0
// 0         0
// 0         sqrt[2]
#define CUP {1,0, 0,0, 266,0, 0,0, 0,237}
#define CAP {1,0,266,0,0, 0,0,0,0,237}

#define BIFMAX 20 //largest fibonacci number we'll ever use
#define MAX 6765  //MAX = fib(BIFMAX)
#define MAXSTRING 1000 // longest string input

int fib[BIFMAX + 1];  // fibonacci numbers
int fibword[MAX];
  // the infinite fibonacci word (A003849 in OEIS).
  // The basis vectors are indexed by Fibbinary numbers
  // (ie. strings of bits that have no adjacent ones, see A003714)
  // and fibword[i] is the last bit of the ith basis vector.

int * prev;        // the matrix as it was at the previous carriage return of ascii input.
int bif_prev_rows; // the number of rows in the "prev" matrix is fib[bif_prev_rows]
int bif_prev_cols;
int * curr;        // the current matrix
int bif_curr_rows;
int bif_curr_cols;
int * next;
int bif_next_rows;
int bif_next_cols;

void die(char* s)
{
  printf("ERROR: %s\n",s);
  free(prev);
  free(curr);
  free(next);
  exit(0);
}

void empty_curr()
// set the curr matrix to represent the empty diagram
{
  curr = malloc(4 * sizeof(int));
  curr[0] = 1;
  curr[1] = 0;
  curr[2] = 0;
  curr[3] = 1;
  bif_curr_rows = 3;
  bif_curr_cols = 3;
}

void initialize()
{
  int i;
  int j;

  //initialize fib[]
  fib[0] = 0;
  fib[1] = 1;
  for(i = 2; i <= BIFMAX; i++)
  {
    fib[i] = fib[i-1] + fib[i-2];
  }

  //initialize prev
  prev = malloc(1 * sizeof(int));
  prev[0] = 1;
  bif_prev_rows = 1;
  bif_prev_cols = 1;

  empty_curr();

  // Initialize fibword to the "fibonacci word" - see the comment when it was declared.
  fibword[0] = 0;
  fibword[1] = 1;
  for(i = 2; i < BIFMAX - 1; i++)
  {
    // Make a copy of the first fib[i] bits shifted by fib[i+1]
    for(j = 0; j < fib[i]; j++)
    {
      if(fib[i + 1] + j >= MAX)
      {
        die("error initializing fibword");
      }
      fibword[fib[i + 1] + j] = fibword[j];
    }
  }
}

void prettyprint()
{
  int i = 0;
  int j = 0;
  int * A = prev;
  int r = fib[bif_prev_rows];
  int c = fib[bif_prev_cols];


  if(r * c > 1000)
  {
    printf("%d by %d is too big to pretty print.\n", r, c);
    return;
  }

  for(i = 0; i < r; i++)
  {
    printf("[ ");
    for(j = 0; j < c; j++)
    {
      printf("%5d ", *A++);
    }
    printf("]\n");
  }
  printf("\n");
}

void multiply()
// prev = prev * curr  (matrix multiplication modulo MOD).
{
  int *X;
  int *Xij;
  int i;
  int j;
  int x;
  int r;
  int c;
  int m;


  if (bif_prev_cols != bif_curr_rows)
  {
    die("tried to multiply matrices with mismatched dimensions.");
  }
  m = fib[bif_prev_cols];

  r = fib[bif_prev_rows];
  c = fib[bif_curr_cols];
  X = malloc(r * c * sizeof(int));
  Xij = X;
  for(i = 0; i < r; i++)
  {
    for(j = 0; j < c; j++)
    {
      *Xij = 0;
      for(x = 0; x < m; x++)
      {
        // add prev[i,x]*curr[x,j]
        *Xij += prev[i*m + x] * curr[x*c + j];
	      *Xij %= MOD;
      }
      Xij++;
    }
  }
  free(prev);
  prev = X;
  bif_prev_cols = bif_curr_cols;
  empty_curr();
}

void tensor()
// curr = curr * next
// where * is the "fibonacci tensor product"
{
  int *X;
  int *Xij;
  int icurr;
  int inext;
  int jcurr;
  int jnext;
  int r;
  int c;

  r = bif_curr_rows + bif_next_rows - 3;
  c = bif_curr_cols + bif_next_cols - 3;

  X = malloc(fib[r] * fib[c] * sizeof(int));
  Xij = X;
  for(icurr = 0; icurr < fib[bif_curr_rows]; icurr++)
  {
    for(inext = (fibword[icurr] ? fib[bif_next_rows - 1] : 0);
        inext < (fibword[icurr] ? fib[bif_next_rows] : fib[bif_next_rows - 1]);
	      inext++)
    {
      for(jcurr = 0; jcurr < fib[bif_curr_cols]; jcurr++)
      {
        for(jnext = (fibword[jcurr] ? fib[bif_next_cols - 1] : 0);
            jnext < (fibword[jcurr] ? fib[bif_next_cols] : fib[bif_next_cols - 1]);
            jnext++)
        {
          *Xij = curr[icurr * fib[bif_curr_cols] + jcurr] * next[inext * fib[bif_next_cols] + jnext];
          *Xij %= MOD;
          Xij++;
        }
      }
    }
  }
  bif_curr_rows = r;
  bif_curr_cols = c;      
  free(curr);
  curr = X;
}

int main()
{
  char tangle[MAXSTRING];
  char line[MAXSTRING];
  int id[9] = {1,0,0, 0,1,0, 0,0,1};
  int cross[25] = CROSS;
  int uncross[25] = UNCROSS;
  int cap[10] = CAP;
  int cup[10] = CUP;
  int i;
  int lines_processed = 0;

  initialize();
  tangle[0]='\0';
  line[0]='\0';
  while(line[0] != '.')
  {
    fgets(line, MAXSTRING, stdin);
    strcat(tangle,line);
  }

  for(i = 0; tangle[i] != '.'; i++)
  {
    switch(tangle[i])
    {
      case ' ' : 
	      break;
      case '|' :
      case '/' :
      case '\\' :
        next = id;
        bif_next_rows = 4;
        bif_next_cols = 4;
        tensor();
        break;
      case '%' :
        next = cross;
        bif_next_rows = 5;
        bif_next_cols = 5;
        tensor();
        break;
      case 'S' :
        next = uncross;
        bif_next_rows = 5;
        bif_next_cols = 5;
        tensor();
        break;
      case 'u' :
	      next = cup;
	      bif_next_rows = 5;
	      bif_next_cols = 3;
	      tensor();
	      break;
      case 'n' :
	      next = cap;
	      bif_next_rows = 3;
	      bif_next_cols = 5;
	      tensor();
	      break;
      case ',' :
      case '\n' :
	      if (lines_processed == 0)
        {
          prev = curr;
	        bif_prev_rows = bif_curr_rows;
	        bif_prev_cols = bif_curr_cols;
	        empty_curr();
	      }
	      else
        {
	        multiply();
	      }
	      lines_processed++;
	      break;
      default :
        printf("What do you mean by %c?\n",tangle[i]);
        die("unknown character");
	      break;
      }
  }
  prettyprint();
  free(prev);
  free(curr);
}
 
