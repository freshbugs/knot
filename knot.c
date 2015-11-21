// C program to input an ascii picture of a tangle
// and output the matrix it would compute in a Fibonacci quantum computer.
//
// We work over the field of integers mod 521, with q=5, so q^5 = -1.
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
#define Q 5

#define QQ ((Q * Q) % MOD)
#define QQQ ((Q * QQ) % MOD)
#define QQQQ ((Q * QQQ) % MOD)
#define PHI ((Q + MOD - QQQQ) % MOD)
#define PHI_INV (PHI - 1)

#define CUP {1,0,  0,0,  PHI_INV,0,  0,0,  0,1}
#define CAP {1,0,1,0,0,  0,0,0,0,PHI}
#define CROSS {MOD-PHI_INV,0,MOD-QQ,0,0,    0,QQQ,0,0,0,    ((MOD-QQ)*PHI_INV)%MOD,0,(QQQQ*PHI_INV)%MOD,0,0,    0,0,0,QQQ,0,    0,0,0,0,MOD-Q}
#define UNCROSS {MOD-PHI_INV,0,QQQ,0,0,    0,MOD-QQ,0,0,0,    (QQQ*PHI_INV)%MOD,0,((MOD-Q)*PHI_INV)%MOD,0,0,    0,0,0,MOD-QQ,0,    0,0,0,0,QQQQ}
#define MERGE {MOD-PHI_INV,0,1,0,0,    0,1,0,0,0,    0,0,0,1,0}
#define SPLIT {MOD-PHI_INV,0,0,    0,1,0,    PHI_INV,0,0,    0,0,1,    0,0,0}


#define BIFMAX 20 //largest fibonacci number we'll ever use
#define MAX 6765  //MAX = fib(BIFMAX)
#define MAXSTRING 1000 // longest string input

int fib[BIFMAX + 1];  // fibonacci numbers
int fibword[MAX];
  // the infinite fibonacci word (A003849 in OEIS).
  // The basis vectors are indexed by Fibbinary numbers
  // (ie. strings of bits that have no adjacent ones, see A003714)
  // and fibword[i] is the last bit of the ith basis vector.

int * next;
int bif_next_rows;
int bif_next_cols;

int * mat[256];
int bif_rows[256];
int bif_cols[256];

void die(char* s)
{
  int i;
  for(i = 0; i < 256; i++)
  {
    if(mat[i] != NULL)
    {
      free(mat[i]);
    }
  }
  if(s[0] != 0)
  {
    printf("ERROR: %s\n",s);
    exit(0);
  }
  exit(1);
}

void empty_curr()
// set the curr matrix mat['@'] to represent the empty diagram
{
  mat['@'] = malloc(4 * sizeof(int));
  mat['@'][0] = 1;
  mat['@'][1] = 0;
  mat['@'][2] = 0;
  mat['@'][3] = 1;
  bif_rows['@'] = 3;
  bif_cols['@'] = 3;
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

  //initialize mat[]
  for(i = 0; i < 256; i++)
  {
    mat[i] = NULL;
  }

  //initialize mat['!']
  mat['!'] = malloc(1 * sizeof(int));
  mat['!'][0] = 1;
  bif_cols['!'] = 1;
  bif_rows['!'] = 1;


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
  int * A = mat['!'];
  int r = fib[bif_rows['!']];
  int c = fib[bif_cols['!']];


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
// mat['!'] *= curr  (matrix multiplication modulo MOD).
{
  int *X;
  int *Xij;
  int i;
  int j;
  int x;
  int r;
  int c;
  int m;


  if (bif_cols['!'] != bif_rows['@'])
  {
    die("tried to multiply matrices with mismatched dimensions.");
  }
  m = fib[bif_cols['!']];

  r = fib[bif_rows['!']];
  c = fib[bif_cols['@']];
  X = malloc(r * c * sizeof(int));
  Xij = X;
  for(i = 0; i < r; i++)
  {
    for(j = 0; j < c; j++)
    {
      *Xij = 0;
      for(x = 0; x < m; x++)
      {
        // add mat['!'][i,x]*mat['@'][x,j]
        *Xij += mat['!'][i*m + x] * mat['@'][x*c + j];
	      *Xij %= MOD;
      }
      Xij++;
    }
  }
  free(mat['!']);
  mat['!'] = X;
  bif_cols['!'] = bif_cols['@'];
  empty_curr();
}

void tensor()
// mat['@'] = mat['@'] * next
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

  r = bif_rows['@'] + bif_next_rows - 3;
  c = bif_cols['@'] + bif_next_cols - 3;

  X = malloc(fib[r] * fib[c] * sizeof(int));
  Xij = X;
  for(icurr = 0; icurr < fib[bif_rows['@']]; icurr++)
  {
    for(inext = (fibword[icurr] ? fib[bif_next_rows - 1] : 0);
        inext < (fibword[icurr] ? fib[bif_next_rows] : fib[bif_next_rows - 1]);
	      inext++)
    {
      for(jcurr = 0; jcurr < fib[bif_cols['@']]; jcurr++)
      {
        for(jnext = (fibword[jcurr] ? fib[bif_next_cols - 1] : 0);
            jnext < (fibword[jcurr] ? fib[bif_next_cols] : fib[bif_next_cols - 1]);
            jnext++)
        {
          *Xij = mat['@'][icurr * fib[bif_cols['@']] + jcurr] * next[inext * fib[bif_next_cols] + jnext];
          *Xij %= MOD;
          Xij++;
        }
      }
    }
  }
  bif_rows['@'] = r;
  bif_cols['@'] = c;      
  free(mat['@']);
  mat['@'] = X;
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
  int merge[15] = MERGE;
  int split[25] = SPLIT;
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
      case 'i' :
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
      case 'h' :
        next = merge;
        bif_next_rows = 4;
        bif_next_cols = 5;
        tensor();
        break;
      case 'y' :
        next = split;
        bif_next_rows = 5;
        bif_next_cols = 4;
        tensor();
        break;
      case ',' :
      case '\n' :
	      if (lines_processed == 0)
        {
          mat['!'] = mat['@'];
	        bif_rows['!'] = bif_rows['@'];
	        bif_cols['!'] = bif_cols['@'];
	        empty_curr();
	      }
	      else
        {
	        multiply();
	      }
	      lines_processed++;
	      break;
      default :
        if(('a' <= tangle[i]) && (tangle[i] <= 'z') && (mat[tangle[i]] != NULL))
        {
          next = mat[tangle[i]];
          bif_rows['!'] = bif_rows[tangle[i]];
          bif_cols['!'] = bif_cols[tangle[i]];
          tensor();
        }
        else
        {
          die("unknown character.");
        }
	      break;
      }
  }
  prettyprint();
  die("");
}
 
