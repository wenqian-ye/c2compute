#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include<iostream>
#include<fstream>
#include<cstring>
//#include<cmath>
#include<vector>

# include <stdio.h>
# include <stdlib.h>

#include<algorithm>
#include<unordered_map>
using namespace std;

ofstream fout_CountPw("kmercount");
ofstream fout_CountPw2("kmerprob");



// In order to give more precise calculation and avoid zero in dominator, 
// we introduce the SCIENTIFIC_NUMBER calculation. 
// We will re-define addtion, multiplication and power in the following sub rountine functions.
// It is contributed by Prof. Minghua Deng, Peking University.

struct SCIENTIFIC_NUMBER
{
	int factor;
	double value;
};


//***change para  : total number of the dataset
int NS = 0;   
int k = 0;  
int d = 0;   // d = 0: take the whole data



// Save the data file directions and names
char* fileDirect[500];
char* fileName[500];


unsigned long power;//parameter, power = 4^(k-1)

//***change para  : the Matrix for save c2, c2star, c2shepp matrix.   
//*******The maximum number of species(datasets) is 200. Users can change the parameters here.
double c2Form[500][500], c2starForm[500][500], c2sheppForm[500][500]; 
    
//receptacle of kmer counts
unordered_map<unsigned long,unsigned long> HashTable;

//***change para : receptacle of kmer counts for all datasets 
unordered_map<unsigned long,unsigned long> HashTableS[500];

//receptacle of Pw(probability of a kmer word) 
unordered_map<unsigned long,SCIENTIFIC_NUMBER > HashPw;

//***change para : receptacle of Pw(probability of a kmer word) for all datasets
unordered_map<unsigned long,SCIENTIFIC_NUMBER > HashPwS[500];

int seqlength=99999;

// seq: save read/genome(A C G T) in each line from file
char seq[99999];

// inverse of seq
char seq_inverse[99999];
//.....................


// Prior[r][m][n][l] records the probability of each kmer word(Pw) in rth dataset. 
// Each kmer word can be decomposed into m As, n Cs, l Gs and k-m-n-l Ts.
// 30 is not related to 30 species. //***change para
// *******The maximum number of k is 29. Users can change this parameter.
SCIENTIFIC_NUMBER Prior[500][30][30][30];  


// records the frequency of A, C, G, T in a dataset
double anum=0, cnum=0, gnum=0, tnum=0, nnum=0;



// Scientific Number calculation funtion:

// TransToReal : trans a SCIENTIFIC NUMBER to a read number
double TransToReal(SCIENTIFIC_NUMBER dSci)
{
    double dReal=0;
    dReal = dSci.value * pow(10,dSci.factor);
    return dReal;

}


// TransToScientific : trans a read number to a SCIENTIFIC NUMBER
SCIENTIFIC_NUMBER TransToScientific(double dReal)
{
	SCIENTIFIC_NUMBER sciTemp;
	int count;

	if( dReal==0.0 )
	{
		sciTemp.factor=0;
		sciTemp.value=0.0;
	}
	else if(dReal>10.0 || dReal<-10.0)
	{
		count=0;
		while(dReal>10.0 || dReal<-10.0)
		{
			dReal /=10.0;
			count++;
		}
		sciTemp.value=dReal;
		sciTemp.factor=count;
	}
	else if( dReal<1.0 && dReal>-1.0)
	{
		count=0;
		while( dReal<1.0 && dReal>-1.0 )
		{
			dReal *=10.0;
			count--;
		}
		sciTemp.value=dReal;
		sciTemp.factor=count;
	}
	else
	{
		sciTemp.value=dReal;
		sciTemp.factor=0;
	}

	return sciTemp;
}


// SciMultiple : Multiplication of two SCIENTIFIC NUMBERS
SCIENTIFIC_NUMBER SciMultiple(SCIENTIFIC_NUMBER left,SCIENTIFIC_NUMBER right)
{
//    cout << "SciMultiple " << endl;
    
	double dTemp;
	SCIENTIFIC_NUMBER sciTemp;
	int count;

	if( left.value==0.0 || right.value==0.0 )
	{
//        cout << "Both 0 " << endl;
        
		sciTemp.value=0.0;
		sciTemp.factor=0;

		return sciTemp;
	}

	// now both left and right element are nonzero
	dTemp=left.value * right.value;
    
//    cout << "left.value " << left.value << endl;
//    cout << "right.value " << right.value << endl;
//    cout << "dTemp " << dTemp << endl;
    
    
	if( dTemp>10.0 || dTemp<-10.0 )
	{
        
//        cout << "10 < dTemp or dTemp < -10 " << endl;
        
		count=0;
		while(dTemp>10.0 || dTemp<-10.0 )
		{
			dTemp /=10.0;
			count++;
		}
		sciTemp.value=dTemp;
		sciTemp.factor=left.factor+right.factor+count;
	}
	else if( dTemp<1.0 && dTemp>-1.0)
	{
//        cout << "dTemp < 1 or dTemp > -1 " << dTemp << endl;
        
		count=0;
		while( dTemp<1.0 && dTemp>-1.0 )
		{
			dTemp *=10.0;
			count--;
		}
        
//        cout << "count " << count << endl;
        
		sciTemp.value=dTemp;
		sciTemp.factor=left.factor+right.factor+count;
        
//        cout << "sciTemp " << sciTemp.value << " " << sciTemp.factor << endl;
	}
	else
	{
        
//        cout << "dTemp normal " << dTemp << endl;
        
		sciTemp.value=dTemp;
		sciTemp.factor=left.factor+right.factor;
	}

	return sciTemp;
}


// SciMultiple : Multiplication between a SCIENTIFIC NUMBERS and a read number
SCIENTIFIC_NUMBER SciMultiple(SCIENTIFIC_NUMBER left,double right)
{
	SCIENTIFIC_NUMBER sciTemp;

	sciTemp=TransToScientific(right);
	sciTemp=SciMultiple(left,sciTemp);

	return sciTemp;
}





// SciAddition : addition of two SCIENTIFIC NUMBERS
SCIENTIFIC_NUMBER SciAddition(SCIENTIFIC_NUMBER left,SCIENTIFIC_NUMBER right)
{
	double dTemp;
	SCIENTIFIC_NUMBER sciTemp;
	int i,count;

	if( left.value==0.0 || right.value==0.0 )
	{
		if( left.value==0.0 )
			return right;
		else
			return left;
	}

	// now the two element are both non zero
	if( left.factor>=right.factor)
	{
		// left element is larger than right element
		dTemp=right.value;
		for(i=0;i<(left.factor-right.factor);i++)
			dTemp /=10.0;
		dTemp +=left.value;
		if( dTemp==0.0 )
		{
			sciTemp.factor=0;
			sciTemp.value=0.0;
			return sciTemp;
		}

		// now dTemp is not zero
		if( dTemp>10.0 || dTemp <-10.0 )
		{
			count=0;
			while(dTemp>10.0 || dTemp<-10.0 )
			{
				dTemp /=10.0;
				count++;
			}
			sciTemp.value=dTemp;
			sciTemp.factor=left.factor+count;
		}
		else if( dTemp<1.0 && dTemp>-1.0 )
		{
			count=0;
			while(dTemp<1.0 && dTemp>-1.0)
			{
				dTemp *=10.0;
				count--;
			}
			sciTemp.value=dTemp;
			sciTemp.factor=left.factor+count;
		}
		else
		{
			sciTemp.value=dTemp;
			sciTemp.factor=left.factor;
		}
		return sciTemp;
	}
	else
	{
		// right element  is larger than left element
		dTemp=left.value;
		for(i=0;i<(right.factor-left.factor);i++)
			dTemp /=10.0;
		dTemp +=right.value;
		if( dTemp==0.0 )
		{
			sciTemp.factor=0;
			sciTemp.value=0.0;
			return sciTemp;
		}

		// now dTemp is not zero
		if( dTemp>10.0 || dTemp <-10.0 )
		{
			count=0;
			while( dTemp>10.0 || dTemp <-10.0 )
			{
				dTemp /=10.0;
				count++;
			}
			sciTemp.value=dTemp;
			sciTemp.factor=right.factor+count;
		}
		else if( dTemp<1.0 && dTemp>-1.0 )
		{
			count=0;
			while(dTemp<1.0 && dTemp>-1.0)
			{
				dTemp *=10.0;
				count--;
			}
			sciTemp.value=dTemp;
			sciTemp.factor=right.factor+count;
		}
		else
		{
			sciTemp.value=dTemp;
			sciTemp.factor=right.factor;
		}
		return sciTemp;
	}
}


// SciAddition : addition between a SCIENTIFIC NUMBERS and a real number
SCIENTIFIC_NUMBER SciAddition(SCIENTIFIC_NUMBER left,double right)
{
	SCIENTIFIC_NUMBER sciTemp;

	sciTemp=TransToScientific(right);
	sciTemp=SciAddition(left,sciTemp);

	return sciTemp;
}

// SciPow : give the power of a scientific number
SCIENTIFIC_NUMBER SciPow(SCIENTIFIC_NUMBER left,double right)
{
	SCIENTIFIC_NUMBER sciTemp;
	double dTemp;
	int iTemp;
/*
	if(left.value==0.0 )
	{
		printf("the base of the power is nagative\n");
		exit(1);
	}
*/
	if(left.value==0.0 )
	{
		sciTemp.factor=0;
		sciTemp.value=0.0;
		return sciTemp;
	}

	dTemp=(log10(fabs(left.value))+left.factor)*right;

	if( dTemp>0.0 )
		iTemp=int(ceil(dTemp));    //ceil(a)是求不小于a的最小整数。floor(a)表示求不大于a的最大整数
	else
		iTemp=int(floor(dTemp));
	sciTemp.factor=iTemp;
	sciTemp.value=pow(10.0,dTemp-iTemp);

	return sciTemp;
}



// The algorithm in SeqKmerCount: to count the frequency of each kmer in a seqence of read/genome
//	k=6
//	ACGTCACGTACGT...
//	ACGTCA index1
//	 CGTCAC index2 = floor(index1/4^5)*4 + 1(C)
//	  GTCACG index3 = floor(index2/4^5)*4 + 2(G)
//	   TCACGT index4 = floor(index3/4^5)*4 + 3(T)
//	    CACGTA index5 = floor(index4/4^5)*4 + 0(A)
//	     ACGTAC index6 = ...
unsigned long SeqKmerCount(char* seq, int k)
{
//..........................
        // int count: The length of char seq 
	int count = 0;     
//..........................
	int i=0,j=0;
	unsigned long index = 0;
	unsigned long total = 0;//total number of the kmers counted in seq
	while(seq[i])
	{
	      //kmer in seq[i,i+1,i+2,...i+k-1] transfered to an index
	      //current index = floor(previous index/4^(k-1))*4 + 0 or 1 or 2 or 3
	      if(seq[i]=='A'|| seq[i] == 'a') {j++; anum++;}
  	      else if(seq[i]=='C'|| seq[i] == 'c') { j++; index++; cnum++;}
              else if(seq[i]=='G'|| seq[i] == 'g') { j++; index+=2; gnum++;}
  	      else if(seq[i]=='T'|| seq[i] == 't') { j++; index+=3; tnum++;}
  	      else { j=0; index=0; nnum++;}//If seq[i] is ambiguous, reset j and index

              if( j == k )
  	      {
  		  HashTable[index]++;
  		  total++;        
  		  index %= power;// current index = floor(previous index/4^(k-1))
  		  j--;//the lengh of seq[i+1,i+2,...i+k-1]
  	      }
  	      index*=4;//current index = floor(previous index/4^(k-1))*4
  	      i++;
	      count++;
        }

//...................

        // Inversed seq from right to left, but not inverse A-T, G-C
	for (int m = 1; m != count + 1; m++)
		seq_inverse[m - 1] = *(seq + count - m);
	seq_inverse[count] = '\0';
	i=0,j=0;
        index = 0;
//...................
        while(seq_inverse[i])
	{
	       // Kmer in seq[i,i+1,i+2,...i+k-1] transfered to an index
	       // Current index = floor(previous index/4^(k-1))*4 + 0 or 1 or 2 or 3
	       if(seq_inverse[i]=='T'|| seq_inverse[i] == 't') {j++;}
  	       else if(seq_inverse[i]=='G'|| seq_inverse[i] == 'g') { j++; index++;}
	       else if(seq_inverse[i]=='C'|| seq_inverse[i] == 'c') { j++; index+=2;}
  	       else if(seq_inverse[i]=='A'|| seq_inverse[i] == 'a') { j++; index+=3;}
  	       else { j=0; index=0; nnum++;}// If seq[i] is ambiguous, reset j and index

	       if( j == k )
  	       {
  		    HashTable[index]++;
  		    index %= power;// Current index = floor(previous index/4^(k-1))
  		    j--;// The lengh of seq[i+1,i+2,...i+k-1]
  	       }
  	       index*=4;// Current index = floor(previous index/4^(k-1))*4
  	       i++;
        }

        return total;      
}





// The algorithm to count Kmer for a dataset

unsigned long DataKmerCount(int DataNum,char* argv_DataNum, int k)
{
    // Initial HashTable
    HashTable.clear();
    HashPw.clear();

    // Initialize total_DataNum. total_DataNum records the total number of kmer in rth dataset
    unsigned long total_DataNum = 0;

    // Initialize anum, cnum, gnum, tnum, recording the frequency of A, T, C, G of rth dataset
    anum=0; cnum=0; gnum=0; tnum=0;

    // Initialize seq[99999999], making 0 at each position   
    for(int ii = 0; ii < seqlength; ii++)
    {
        seq[ii]=0;
    }
        
    // Parameters: k and power = 4^(k-1)
    power = 1; for( int i = 0; i < k-1; i++) power *= 4;   

	// Open Data from file argv[r*2-1]
	char combinefile[1000];
	strcpy_s(combinefile, fileName[DataNum]);
	strcat_s(combinefile, "_combine");
	cout << "    count kmer from the file: " << combinefile << endl;

	ifstream fin(combinefile);   //EDIT fin(argv[1])


								 //  unsigned long copynumber;
								 //  cout << "begin scan" << endl;

    
    // Begin to scan lines in dataset from top to the bottom

	//这里做了修改本来是全读的
  while(fin.getline(seq,seqlength)) 
    {
       total_DataNum += SeqKmerCount(seq,k);
//    cout << "seq" << endl;
    }
	//现在只读前一部分
	//for (int l_num=1;l_num<=3500;l_num++)
	//{
		
			

		//fin.getline(seq,seqlength);
		//total_DataNum += SeqKmerCount(seq,k);
	//}
    fin.close();


     // Sort the kmer
     vector<unsigned long> temp;
     for( unordered_map<unsigned long,unsigned long>::iterator i = HashTable.begin(); i!= HashTable.end(); i++) temp.push_back(i->first);
     sort(temp.begin(), temp.end());


     unsigned long key;
     for(vector<unsigned long>::iterator j = temp.begin(); j!=temp.end(); j++)
     {
         key = *j;
	 HashTableS[DataNum][key]=HashTable[key];
//	 fout2.write((char*)&key,sizeof(unsigned long));
//	 fout2.write((char*)&HashTableS[DataNum][key],sizeof(unsigned long));
//	 fout2 << "Key" << key;
//	 fout2 << key << "," << HashTable[key] << endl;
//	 fout2 << "Count" << HashTableS[DataNum][key] << endl;
     }
//	fout2.close();


     return total_DataNum;

}



// Compute Pw for each kmer word
void ComputePw(int DataNum, int k, char* argg_DataNumK, char* argg_DataNum, int d)
{
     
     //Compute probability(frequency) of A, C, G, T: pa, pc, pg, pt

     // Initialize pa, pc, pg, pt 
     SCIENTIFIC_NUMBER pa, pc, pg, pt;
     pa.value = 0; pc.value = 0; pg.value = 0; pt.value = 0;
     pa.factor = 0; pc.factor = 0; pg.factor = 0; pt.factor = 0;

     // Total count of A, C, G, T
     double totalACGT = 0;
     totalACGT = anum + cnum + gnum + tnum;

     // Frequency of A, C, G, T
     pa=TransToScientific(anum/totalACGT);
     pc=TransToScientific(cnum/totalACGT);
     pg=TransToScientific(gnum/totalACGT);
     pt=TransToScientific(tnum/totalACGT);
//   cout<<pa<<pc<<pg<<pt<<endl;

 
     // Output probability of A, C, G, T
     char dstr[5]; sprintf_s(dstr, "%d", d); 
     char outputfile_pa[1000];
     strcpy_s(outputfile_pa,argg_DataNum);
     strcat_s(outputfile_pa,"_k");
     strcat_s(outputfile_pa,argg_DataNumK);
     strcat_s(outputfile_pa,"_d"); 
     strcat_s(outputfile_pa,dstr);
     strcat_s(outputfile_pa,"_pa"); // strcpy(outputfile,argv[1]) strcat(outputfile,argv[2])
     ofstream fout_pa(outputfile_pa);

     // Convert SCITIFIC_NUMBER to real number
     long double pa_real = TransToReal(pa);
     long double pc_real = TransToReal(pc);
     long double pg_real = TransToReal(pg);
     long double pt_real = TransToReal(pt);
     fout_pa << "pa" << pa_real << endl;
     fout_pa << "pc" << pc_real << endl;
     fout_pa << "pg" << pg_real << endl;
     fout_pa << "pt" << pt_real << endl << endl;
     fout_pa.close();


     // Computer Pw(probability) for all kmer word

     // First compute A[] C[] G[] T[]: A[pa^0, pa^1, pa^2, ... , pa^k] 
	 SCIENTIFIC_NUMBER *A = new SCIENTIFIC_NUMBER[k+1],*C = new SCIENTIFIC_NUMBER[k+1], *G = new SCIENTIFIC_NUMBER[k+1], *T = new SCIENTIFIC_NUMBER[k+1];
     for (int p = 0; p < k + 1; p++)
     {
//         cout << "p" << p << " k" << k << endl;
         
         A[p]=SciPow(pa,p);
//         cout<< "A[p]" << A[p].value << A[p].factor << endl;
         C[p]=SciPow(pc,p);
//         cout<< "C[p]" << C[p].value << C[p].factor << endl;
         G[p]=SciPow(pg,p);
//         cout<< "G[p]" << G[p].value << G[p].factor << endl;
         T[p]=SciPow(pt,p);
//         cout<< "T[p]" << T[p].value << T[p].factor << endl;
     }
//     cout<<endl;


     // Second calculate probability for each kmer word: 
     // each kmer word can be decomposed into n_1 As, n_2 Cs, n_3 Gs, and n_4 Ts, such that n_1 + n_2 + n_3 + n_4 = k
     int m, n, l;
     for( m = 0; m <k+1; m++)
     {
      //   Prior[m][n][l] = A[m];
         for( n = 0; n < k+1 - m; n++)
         {
          //   Prior[m][n][l]*=C[n];
             for( l = 0; l < k+1 -m - n; l++)
             {
                    Prior[DataNum][m][n][l] = SciMultiple(SciMultiple(SciMultiple(A[m], C[n]), G[l]), T[k-m-n-l]);
//                    cout << "m" << m << "n" << n << "l" << l << " Prior " <<  Prior[DataNum][m][n][l].value << " " << Prior[DataNum][m][n][l].factor <<endl;
//                    cout << "Prior-020 " << Prior[DataNum][0][2][0].value << " " << Prior[DataNum][0][2][0].factor <<endl;
             }
         }
     }
//   cout<<endl;
	 delete A,C,G,T;


}


// Print kmer count hashtable and its corresponding Pw
void PrintKmerCountPw(int k, int DataNum, char* argg_DataNumK, char* argg_DataNum, int d)
{
     
     // Output kmer count-pw files
     char dstr[5]; sprintf_s(dstr, "%d", d); 
     char outputfile_CountPw[1000]; 
     strcpy_s
	(outputfile_CountPw,argg_DataNum); 
     strcat_s(outputfile_CountPw,"_k");
     strcat_s(outputfile_CountPw,argg_DataNumK);
     strcat_s(outputfile_CountPw,"_d"); // strcpy(outputfile,argv[1]) strcat(outputfile,argv[2])
     strcat_s(outputfile_CountPw,dstr);
     strcat_s(outputfile_CountPw,"_wordcount_pw"); // strcpy(outputfile,argv[1]) strcat(outputfile,argv[2])
     //ofstream fout_CountPw(outputfile_CountPw);


     for(unsigned long kk = 0; kk < pow(4,k); kk++)
     {

         // For each kmer, convert base-10 numeral system back to base-4 numeral system
		 int a; int b[30]; memset(b, 0, sizeof(b));
         int t=k-1;
         a = kk; 
         while(a>=3) {b[t]=a%4;a/=4;t--;}
         b[t]=a;   //4 jinzhi

//       for(int i=0;i<k;i++)
//       {
//           cout<<b[i]<<endl;
//       }

         // Count number of As Cs Gs Ts in a kmer word
         int count_ACGT[4]={0, 0, 0, 0};
         for (int mm = 0; mm < k; mm++)
         {
             count_ACGT[b[mm]]++;
         }
//         cout << "count " << count_ACGT[0] << count_ACGT[1] << count_ACGT[2] << count_ACGT[3] << endl;


         // Search specific kmer word probability
         long double pwREAL = TransToReal(Prior[DataNum][count_ACGT[0]][count_ACGT[1]][count_ACGT[2]])+TransToReal(Prior[DataNum][count_ACGT[3]][count_ACGT[2]][count_ACGT[1]]);

      
         // Output count-pw files              
         fout_CountPw <<HashTableS[DataNum][kk] <<"\t" ;
		 fout_CountPw2 <<pwREAL <<"\t" ;
//         cout << kk << "," << HashTableS[DataNum][kk] << "," << pwREAL << endl;
         
			 
              
      }
	 fout_CountPw<<endl;
	 fout_CountPw2<<endl;

			//fout_CountPw.close();
}



void D2C2compute(int k, unsigned long total[])
{

//    cout << "k1" << k << endl;

    int knew = k;

    for (int pp = 1; pp <= NS; pp++)
    {      
         c2Form[pp][pp] = 1; 
         c2starForm[pp][pp] = 1; 
         c2sheppForm[pp][pp] = 1; 
      
        for (int qq = pp + 1; qq <= NS; qq++)
        {
            
            //cout << "Compute D2C2 value for Dataset No." << pp << " and No." << qq << "." << endl;
            //cout << "k2 " << k << " knew " << knew << endl;
            
            SCIENTIFIC_NUMBER wordmean1, wordmean2;
            wordmean1.value = 0.0; wordmean1.factor = 0;
            wordmean2.value = 0.0; wordmean2.factor = 0;

            SCIENTIFIC_NUMBER tilde1, tilde2;
            tilde1.value = 0.0; tilde1.factor = 0;
            tilde2.value = 0.0; tilde2.factor = 0;

            SCIENTIFIC_NUMBER tildeprod;
            tildeprod.value = 0.0; tildeprod.factor = 0;

            SCIENTIFIC_NUMBER vc;
            vc.value = 0.0; vc.factor = 0;

            SCIENTIFIC_NUMBER x2, y2;
            x2.value = 0.0; x2.factor = 0;
            y2.value = 0.0; y2.factor = 0;

            SCIENTIFIC_NUMBER x2tail, y2tail;
            x2tail.value = 0.0; x2tail.factor = 0;
            y2tail.value = 0.0; y2tail.factor = 0;

            SCIENTIFIC_NUMBER x2s, y2s;
            x2s.value = 0.0; x2s.factor = 0;
            y2s.value = 0.0; y2s.factor = 0;

            SCIENTIFIC_NUMBER d2, d2star, d2shepp,manhattan,euclid;
            d2.value = 0.0; d2.factor = 0;
            d2star.value = 0.0; d2star.factor = 0;
            d2shepp.value = 0.0; d2shepp.factor = 0;
			manhattan.value=0.0;manhattan.factor=0;
			euclid.value=0.0;euclid.factor=0;

            SCIENTIFIC_NUMBER c2_deno, c2star_deno, c2shepp_deno;
            c2_deno.value = 0.0; c2_deno.factor = 0;
            c2star_deno.value = 0.0; c2star_deno.factor = 0;
            c2shepp_deno.value = 0.0; c2shepp_deno.factor = 0;

            SCIENTIFIC_NUMBER c2, c2star, c2shepp;
            c2.value = 0.0; c2.factor = 0;
            c2star.value = 0.0; c2star.factor = 0;
            c2shepp.value = 0.0; c2shepp.factor = 0;


/*
            //output d2c2
            char outputfile8[1000];
            strcpy(outputfile8,argv[pp*2-1]);
            strcat(outputfile8,"_");
            strcat(outputfile8,argv[qq*2-1]);
            strcat(outputfile8,"_k");
            strcat(outputfile8,argv[pp*2]);
            strcat(outputfile8,"_d2c2");
//          cout << outputfile8 <<endl;
	    ofstream fout(outputfile8);
*/
           



			SCIENTIFIC_NUMBER  normalize1,normalize2,inv_nor1,inv_nor2;
			normalize1.value=0.0; normalize1.factor=0;
			normalize2.value=0.0; normalize2.factor=0;

			for(unsigned long kk = 0; kk < pow(4,knew); kk++)
            {

                
				int a; int b[30]; memset(b, 0, sizeof(b));
                a = kk; 
                int t=knew-1;
                while(a>=3) {b[t]=a%4;a/=4;t--;}
                b[t]=a;   //4 jinzhi



                //count A T C Gs in kmer word
                int count[4]={0, 0, 0, 0};
                for (int mm=0; mm<knew; mm++)
                {
                     count[b[mm]]++;
                }

                SCIENTIFIC_NUMBER Kmercount1, Kmercount2;
                Kmercount1 = TransToScientific(HashTableS[pp][kk]);
                Kmercount2 = TransToScientific(HashTableS[qq][kk]);
              
				normalize1=SciAddition(Kmercount1,normalize1);
				normalize2=SciAddition(Kmercount2,normalize2);
				

        }

			inv_nor1.value=1/normalize1.value; inv_nor1.factor=-normalize1.factor;
			inv_nor2.value=1/normalize2.value; inv_nor2.factor=-normalize2.factor;

            
            // Compute Pw for each Kmer from Prior[][][][]



            for(unsigned long kk = 0; kk < pow(4,knew); kk++)
            {
//                cout << "kk " << kk << " knew " << knew << endl;
                
                
				int a; int b[30]; memset(b, 0, sizeof(b));
                a = kk; 
                int t=knew-1;
                while(a>=3) {b[t]=a%4;a/=4;t--;}
                b[t]=a;   //4 jinzhi


                //count A T C Gs in kmer word
                int count[4]={0, 0, 0, 0};
                for (int mm=0; mm<knew; mm++)
                {
                     count[b[mm]]++;
                }              

                SCIENTIFIC_NUMBER Kmercount1, Kmercount2;
                Kmercount1 = TransToScientific(HashTableS[pp][kk]);
                Kmercount2 = TransToScientific(HashTableS[qq][kk]);



                SCIENTIFIC_NUMBER xhat,yhat,negyhat,wordman,abswordman,wordeuc;
				xhat=SciMultiple(Kmercount1,inv_nor1);
				yhat=SciMultiple(Kmercount2,inv_nor2);
				negyhat.value=-yhat.value;
				negyhat.factor=yhat.factor;
				wordman=SciAddition(xhat,negyhat);
				abswordman.value=abs(wordman.value);
				abswordman.factor=wordman.factor;
				wordeuc=SciPow(SciAddition(xhat,negyhat),2);

             
				manhattan=SciAddition(manhattan, abswordman);
				euclid=SciAddition(euclid,wordeuc);

              

            }
			euclid=SciPow(euclid,0.5);


            double manReal,eucReal;

			manReal=TransToReal(manhattan);
			eucReal=TransToReal(euclid);


            c2Form[pp][pp] = 1; c2Form[qq][pp] = manReal;
            c2starForm[pp][pp] = 1; c2starForm[qq][pp] = eucReal;
            c2sheppForm[pp][pp] = 1; c2sheppForm[qq][pp] = 1;


        }

    }


}




void PrintD2C2(int d, char* argv_DataNumK)
{

     char dstr[5]; sprintf_s(dstr, "%d", d); 
       
     // output csv file
     char outputfile_c2[1000];
     strcpy_s(outputfile_c2,"k");
     strcat_s(outputfile_c2,argv_DataNumK);
     strcat_s(outputfile_c2,"-d");
     strcat_s(outputfile_c2,dstr);
     strcat_s(outputfile_c2,"-c2Form");
//     cout << outputfile_c2 << endl;
     ofstream fout_c2(outputfile_c2);

     for(int t1 = 1; t1 <= NS; t1++)
     {
          for(int t2 = 1; t2 <= NS; t2++)
          {
               if(t2 < t1) { fout_c2 << c2Form[t1][t2] << "\t" ; }
               else if(t2 >=t1 && t2 != NS ) { fout_c2 << "\t"; }
               else {fout_c2 << "\n"; }  
           }
     }
       
     fout_c2.close();




      char outputfile_c2star[1000];
      strcpy_s(outputfile_c2star,"k");
      strcat_s(outputfile_c2star,argv_DataNumK);
      strcat_s(outputfile_c2star,"-d");
      strcat_s(outputfile_c2star,dstr);
      strcat_s(outputfile_c2star,"-c2starForm");
      ofstream fout_c2star(outputfile_c2star);
           
      for(int t1 = 1; t1 <= NS; t1++)
      {
           for(int t2 = 1; t2 <= NS; t2++)
           {
                if(t2 < t1) { fout_c2star << c2starForm[t1][t2] << "\t" ; }
                if(t2 >=t1 && t2 != NS ) { fout_c2star << "\t"; }
                if(t2==NS){fout_c2star << "\n"; }                     
           }
      }
       
      fout_c2star.close();



     char outputfile_c2shepp[1000];
     strcpy_s(outputfile_c2shepp,"k");
     strcat_s(outputfile_c2shepp,argv_DataNumK);
     strcat_s(outputfile_c2shepp,"-d");
     strcat_s(outputfile_c2shepp,dstr);
     strcat_s(outputfile_c2shepp,"-c2sheppForm");
     ofstream fout_c2shepp(outputfile_c2shepp);


     for(int t1 = 1; t1 <= NS; t1++)
     {
          for(int t2 = 1; t2 <= NS; t2++)
          {
                    
               if(t2 <= t1) {  fout_c2shepp <<c2sheppForm[t1][t2] << "\t" ; }
               if(t2 >t1 && t2 != NS ) { fout_c2shepp << "\t"; }
               if(t2==NS){fout_c2shepp << "\n"; }                     
          }
     }
       
     fout_c2shepp.close();


}





// Extract the raw data(.fasta/fastq) file name and its direction. 
int AllocateFile(char* argv_DataFileName){

        char seq[500][999];
        int seqlength=999; 
        char delims[] = " ";
        char *dir_name = NULL;
      
 
        ifstream infile(argv_DataFileName); 
        int lineNum = 1;

        while(infile.getline(seq[lineNum],seqlength))
        {	

			char *buffffff= new char[1000];
            dir_name = strtok_s( seq[lineNum], delims, &buffffff);

            int i = 0;
            while( dir_name != NULL ){
                
                if(i == 0){
                    
                    fileDirect[lineNum] = dir_name; 
                         
                    dir_name = strtok_s( NULL, delims, &buffffff);
                    i ++;
                    continue; 
                                      
                }
                
                if(i == 1){
                        
                    fileName[lineNum] = dir_name; 
                    break;    
                }
            }

            
            lineNum ++;
            

   
       }
    
       

       return (lineNum-1);


}
    
    

// Data Pre-process: Transform .fasta(/.fastq) to read sequence only data.
void dataPreProcess(int DatasetNum, char *fastType){


    char outputfile[1000];
    strcpy_s(outputfile,fileName[DatasetNum]);
    strcat_s(outputfile,"_combine");
    ofstream fout(outputfile); 
 

    char seq[99999];
    int seqlength=99999;
    ifstream infile(fileDirect[DatasetNum]);
    

    if(fastType == "Q" || fastType == "q"){
         
         cout << seq << endl;
        
         int flagA = 0;
         while(infile.getline(seq,seqlength))
         {
             
             
    	       if(seq[0]=='>')
    	       {
    	             if(flagA == 1){fout << endl;}
                     continue;
    	       }
    	       else
    	       {
                     flagA = 1;
                     fout << seq;
               }
         }
    }
    else{
         
         int flagQ = 0;
        
         while(infile.getline(seq,seqlength))
         {
             
    	       if(seq[0]=='A' || seq[0]=='a' || seq[0]=='C' || seq[0]=='c' || seq[0]=='G' || seq[0]=='g' || seq[0]=='T' || seq[0]=='t' || seq[0]=='N' || seq[0]=='n')
    	       {
                     flagQ = 1;
    	             fout<<seq;
    	       }
    	       else
    	       {
                     if(flagQ == 1){fout << endl;}
                     flagQ = 0;
                     continue;
               }
         }  


   }


}









//KmerCount.out [sample data file] [k]
int main(int argc, char *argv[])   //EDIT main(int argc, char *argv[])
{

    // Extract the User-specific parameters.
	argv[1] = "6";
	argv[2] = "F:/XinBai/livercancerproject/Data/testfree";
	argv[3] = "A";

    k = atoi(argv[1]);

    NS = AllocateFile(argv[2]);
    
    char *dataType = argv[3];
    

    HashTable.clear();
    HashPw.clear();
	unsigned long *total = new unsigned long[NS+1];     //*****change parameters

//    cout << fileName[1] << fileName[2] << fileName[3] << endl;

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Data Pre-process: Eliminate the annotation lines in raw dataset .fastq/.fasta. Only leave the read sequences
    cout << "0. Data Pre-process: Eliminate the annotation lines in raw dataset .fastq/.fasta. Only leave the read sequences..." << endl;
    
    for (int DataNum = 1; DataNum <= NS; DataNum++)     //////
    {
          
        cout << "Data Pre-process Dataset No." << DataNum << "." << endl;
        
        dataPreProcess(DataNum, dataType);
 
    }





    // Count frequency for each Kmer word and compute Pw for each Kmer word in each dataset
    cout << "1. Count frequency for each Kmer word and compute Pw for each Kmer word in each dataset..." << endl;
    
    for (int DataNum = 1; DataNum <= NS; DataNum++)     //////
    {
          
        cout << "Scanning Dataset No." << DataNum << "." << endl;
        
         // Count kmer in rth dataset
         total[DataNum] = DataKmerCount(DataNum, fileDirect[DataNum], k);

         // Compute Probability of each kmer word(Pw)
         ComputePw(DataNum, k, argv[1], fileName[DataNum], d);
 
    }




    // Output key, hashtable and pw for each dataset
    
    cout << "2. Print out key, hashtable and pw for each dataset..." << endl;
 
    for (int DataNum = 1; DataNum <= NS; DataNum++)  
    {
        
         cout << "Printing out Kmer count and its Pw for Dataset No." << DataNum << "." << endl;
		 

         // Print kmer count hashtable and its corresponding Pw
         PrintKmerCountPw(k, DataNum, argv[1],fileName[DataNum], d);

    }
	fout_CountPw.close();
	fout_CountPw2.close();


    // Compute D2C2 for any pair of datasets from their Kmer-Pw information 
    
    cout << "3. Compute D2C2 for any pair of datasets from their Kmer-Pw information..." << endl;

//    cout << "k0" << k << endl;
    
    D2C2compute(k, total);


    
    // Output D2C2 result matrix: c2Form, c2starForm, c2sheppForm
    
    cout << "4. Print D2C2 result matrix: c2Form, c2starForm, c2sheppForm..." << endl;

    PrintD2C2(d, argv[1]);



	delete total;
	system("PAUSE");
    return 0;
}






/////////////////////////////////////END/////////////////////////////////////////////////////

