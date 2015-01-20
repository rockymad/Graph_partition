#include<iostream>
#include<vector>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<stdexcept>
#include <cmath>
using namespace std;
unsigned short VERBOSITY=1;



#define myabs(a) (((a) > 0) ? (a):(-(a)))

unsigned int size();
void setDoEigenvectors(bool doit);
void setDiagonalizationInterval(unsigned int interval) ;
  void setMaxIterations(unsigned int iters);
  void setEpsilon(float eps);
  vector<vector<float> > getEigenvectors();
    vector<float> getEigenvalues();
    std::vector<float> readVectorFromFile(std::string);
    std::vector<int> readVectorFromFilec(std::string);
    void writeVector(std::vector<float> vec,std::string name="");
    void writeVector(float* vec,unsigned int length,std::string name="");
    void writeNonZeroEntries(float* vec, unsigned int length,std::string name="");   
    void solver(std::vector<float>, std::vector<float>);
    void Doo();
    void NRjacobi(float **, int, int , float **, int *);
    void calculateEigenvectors(std::vector<float> ,std::vector<std::vector<float> > );
    void sortEigenvalues();
	void jacobi_method ( float **, float **, int);
	float maxoffdiag ( float ** A, int * k, int * l, int n );
		void rotate ( float ** A, float ** R, int k, int l, int n );
    float **matrix_;
    bool delete_matrix_;
    



unsigned int n_;
float matrix[20][20];
void loadFromFiles(string fileprefix);
void readarr(string filename);
void write();
void cpu_multiply(float* in, float *out);

void loadTestMatrix();
unsigned int matrix_size_;
unsigned int max_cnt_;
int *dev_col_;
float *dev_ele_;
std::vector<int> cnt_;
std::vector<int> col_;
std::vector<float> ele_;
extern unsigned short VERBOSITY;
void Initialize(std::vector<float> init_vec);
void DoStep();
void free_vecs();
void calculateEigenvectors(std::vector<float> init_vec,std::vector<std::vector<float> > reduced_eigenvectors);
unsigned int grid_;
unsigned int block_;
float *prev_;
float *act_;
float *next_;
 static const unsigned int DEFAULT_MAX_ITERATIONS=1000000;
  static const float DEFAULT_EPSILON=1E-3;
  void Do();
  void Initialize(std::vector<float> init_vec);
  float residual(float val, std::vector<float> vec);
void calculateEigenvectors(std::vector<float> init_vec,std::vector<std::vector<float> > reduced_eigenvectors);
  unsigned int max_iterations_;
  unsigned int number_eigenvalues_;
  unsigned int step_;
  vector<float> k_;
  vector<float> e_;
  bool do_eigenvectors_;
  unsigned int diagonalization_interval_;
  vector<float> eigenvalues_;
  vector<vector<float> > eigenvectors_;
  float epsilon_;



  








int main(int argc, char** argv)
{
  std::string add_graph="example";
  
  dev_col_=NULL;
  dev_ele_=NULL;
  
  loadFromFiles(add_graph);
  readarr(add_graph);
  cout<<"Matrix size: "<<size()<<endl;
    cout<<"CPU"<<endl;
      max_iterations_=DEFAULT_MAX_ITERATIONS;
      do_eigenvectors_=false;
      number_eigenvalues_=2;
      diagonalization_interval_=2;
      epsilon_=DEFAULT_EPSILON;
      prev_=NULL;
      act_=NULL;
      next_=NULL;
    setMaxIterations(500);
    setEpsilon(1E-3);
    setDiagonalizationInterval(20);
    setDoEigenvectors(true);
    Do();
    std::vector<float> cpu_vals=getEigenvalues();
    for(unsigned int i=0;i<cpu_vals.size();i++)
    cout<<std::setprecision(10)<<cpu_vals[i]<<" "<<endl;
    
    std::vector<std::vector<float> > cpu_vecs=getEigenvectors();
 
}


void readarr(string filename)
{
	ifstream ifs(filename.c_str(),std::ios::in);	
	  
	  if (!ifs.is_open())
	      {
		  cout<<"File "+filename+" could not be opened."<<endl;
	      }
	  string line;	
	  int i=0;
	  while(!ifs.eof())
	    { 
		  int j=0;
	      getline(ifs,line);
	       if(line=="")
		   break;
	      stringstream sline(line);
	      string str1;
	      while(sline>>str1)
	      {
	      float val;
	      stringstream strr(str1);
	      while(strr>>val)
	      {
	    	  matrix[i][j]=val;
	    	 // cout<<val<<endl;
	    	  j++;
	      }
	      }
	      i++;
	      //cout<<val;

	    }
	  ifs.close();
}


void loadFromFiles(string fileprefix)
{

  cnt_=readVectorFromFilec(fileprefix+"_CNT.dat");
  col_=readVectorFromFilec(fileprefix+"_COL.dat");
  ele_=readVectorFromFile(fileprefix+"_ELE.dat");
  matrix_size_=cnt_.size();
  cout<<"ele:"<<ele_.size()<<" cnt"<<cnt_.size()<<" col"<<col_.size()<<endl;
  if(col_.size()!=ele_.size())
    throw runtime_error("loadMatrix number_elements_");
}


void write()
{
  cout<<"Matrix size: "<<matrix_size_<<endl;
  cout<<"Max cnt: "<<max_cnt_<<endl;
  unsigned int data_len=matrix_size_*max_cnt_;
  int* col=new int[data_len];
  for(unsigned int i=0;i<data_len;i++)
    cout<<col[i]<<endl;
  float* ele=new float[data_len];
  for(unsigned int i=0;i<data_len;i++)
    cout<<ele[i]<<endl;

  delete[] col;
  delete[] ele;
  
   
}

 void cpu_multiply(float* in, float *out) 
 {
	 
	 
   /*unsigned  int ele_cnt=0;
   #pragma omp parallel for shared(matrix_size_)
   omp_set_num_threads(NUM_THREADS);
   for(unsigned int row=0;row<matrix_size_;row++)
     {
       out[row]=0.0;
       for(unsigned int entry=0;entry<cnt_[row];entry++)
	 {
	   out[row]+=ele_[ele_cnt]*in[col_[ele_cnt]];
	   ele_cnt++;
	 }
     }
*/

	   for(unsigned int row=0;row<size();row++)
	     {
		   out[row]=0.0;
	     for(unsigned int entry=0;entry<size();entry++)
		 {
		   out[row]+=matrix[row][entry]*in[entry];
		   
		 }

	     }
   
 }
void loadTestMatrix()
{
  unsigned int number_elements=5;
  unsigned int matrix_size=3;
  cnt_.resize(matrix_size);
  col_.resize(number_elements);
  ele_.resize(number_elements);
  cnt_[0]=1;  cnt_[1]=2;  cnt_[2]=2;
  col_[0]=1;  col_[1]=0;  col_[2]=2;  col_[3]=1;  col_[4]=2;
  ele_[0]=2.0;ele_[1]=2.0;ele_[2]=5.0;ele_[3]=5.0;ele_[4]=4.0;
 
}
void Do()
{
  //time_t time=clock();
  std::vector<float> initial_vec(size(),0.0);
  initial_vec[0]=1.0;
  e_.clear();
  k_.clear();
  Initialize(initial_vec);
  std::vector<std::vector<float> > reduced_eigenvectors;
  cout<<"Lanczos computations ->"<<endl;
  time_t time1=clock();
  for(unsigned int step=0;step<500;step++)
    {  
      DoStep();
      if(step % diagonalization_interval_==0 && step>0)
  	{
      	 
      	  
  	  //TridiagonalSolver solver(e_,k_);
  	  
  	  solver(e_,k_);
  	  setMaxIterations(1000);
  	  Doo();
  	  sortEigenvalues();
  	vector<float> vals=getEigenvalues();
  		  
  		  cout<<"Eigenvalues in step "<<step<<": ";
  		  //for(unsigned int i=0;i<vals.size();i++)
  		  //  cout<<vals[i]<<" ";
  		  cout<<endl;
  		    
  	 	  if(eigenvalues_.size()>=number_eigenvalues_ )
  		    {
  		      bool converged=false;
  		      for(unsigned int i=0;i<number_eigenvalues_;i++)
  		      {
  		    	//  cout<<"epsilon"<<epsilon_<<endl;
  		    	  //cout<<"fabs(vals[i]-eigenvalues_[i] ="<<vals[i]<<"--"<<eigenvalues_[i]<<"===="<<myabs(vals[i])-myabs(eigenvalues_[i])<<endl;
  			if(step>=200)
  			  converged=true;
  		      }
  		      if(converged)
  			{
  			  cout<<"converged after "<<step<<" iterations."<<endl;
  			  eigenvalues_=vals;
  			  reduced_eigenvectors=getEigenvectors();
  			  eigenvalues_.resize(number_eigenvalues_);
  			  reduced_eigenvectors.resize(number_eigenvalues_);
  			  break;
  			}
  		    }
  		  eigenvalues_=vals;
  	}
      if(step==max_iterations_-1)
	{
	  throw std::runtime_error("Max iterations reached, Lanczos not converged.");
	}
    }
  cout<<"primary iterations completed in: "<<(float)(clock()-time1)/CLOCKS_PER_SEC<<"s."<<endl;
   //eigenvectors
   if(do_eigenvectors_)
     {
       cout<<"calculating eigenvectors..."<<endl;
       time_t time2=clock();
       
      
       
       calculateEigenvectors(initial_vec,reduced_eigenvectors);
     //  cout<<"eigenvectors calculated in: "<<(float)(clock()-time2)/CLOCKS_PER_SEC<<"s."<<endl;
     }

   //cout<<" Lanczos iterations completed in: "<<(float)(clock()-time)/CLOCKS_PER_SEC<<"s."<<endl;
 
}
void Initialize(std::vector<float> init_vec) {
	cout << "Initialization..." << endl;
	time_t time = clock();
	
	if (init_vec.size() != size())
		throw runtime_error("GPULanczos::Initialize wrong start vector size");
	
	next_ = new float[size()];  //10
	act_ = new float[size()];
	prev_ = new float[size()];

	for (unsigned int i = 0; i < size(); i++) {
		next_[i] = init_vec[i];
		act_[i] = 0.0;
		prev_[i] = 0.0;
	}

	step_ = 0;
	cout << "completed in " << (float) (clock() - time) / CLOCKS_PER_SEC
			<< "s." << endl;
}
void DoStep()
{
	int check;
  if(VERBOSITY>2)
    cout<<"iterating step "<<step_<<endl;
  if (step_>0)
    {
      for(unsigned int i=0;i<size();i++)
	next_[i]-= e_[step_-1]*act_[i];

      if (step_>1)
	for(unsigned int i=0;i<size();i++)
	  next_[i]-= k_[step_-2]*prev_[i];


      if(k_.size()<=step_-1)
	{
	 float norm=0.0;
	  for(unsigned int i=0;i<size();i++)
	  {
	    norm+=next_[i]*next_[i]; 
	 
	  }
	  k_.push_back(sqrt(norm));
	}
	if(VERBOSITY>1)
	  cout<<"k:"<<k_[step_-1]<<endl;
      //FIX check if k too small
     
      for(unsigned int i=0;i<size();i++)
	next_[i]/=k_[step_-1];

   }
 
  float *temp=prev_;
  prev_=act_;
  act_=next_;
  next_=temp;
  
/*  cout<<"The act_ look like ";
  for(int i=0;i<10;i++)
	  {
	    cout<<act_[i];
	  }
  cout<<endl<<"The next_ look like ";
  for(int i=0;i<10;i++)
  	  {
  	    cout<<next_[i];
  	  }
  cout<<endl;*/
  cpu_multiply(act_,next_);
/* cout<<"The act_ look like ";
  for(int i=0;i<10;i++)
	  {
	    cout<<act_[i];
	  }
  cout<<endl<<"The next_ look like ";
  for(int i=0;i<10;i++)
  	  {
  	    cout<<next_[i];
  	  }
  cout<<endl;*/
  //cin>>check;
cout<<"The size of e_ ="<<e_.size()<<endl;
  if(e_.size()<=step_)
    {
      float prod=0.0;
      for(unsigned int i=0;i<size();i++)
      {
	prod+=act_[i]*next_[i];
      }
      cout<<"E_val "<<prod<<endl;
      e_.push_back(prod);//0
    }
  if(VERBOSITY>1)
    cout<<"e:"<<e_[step_]<<endl;
  step_++;
  //cout<<"step finished"<<endl;
  //writeDeviceVector(dev_next_,size(),"next");
}


/////IO////

std::vector<float> readVectorFromFile(std::string filename)
{		    
  std::ifstream ifs(filename.c_str(),std::ios::in);	
  
  if (!ifs.is_open())
      throw std::runtime_error("File "+filename+" could not be opened.");

  std::string line;			   
  std::vector<float> vec(0);
  while(!ifs.eof())
    { 
      getline(ifs,line);
       if(line=="")
	break;
      std::stringstream sline(line);
      float val;
      sline>>val;
      vec.push_back(val);
    }
  ifs.close();

  cout<<"Vector "<<filename<<" read, "<<vec.size()<<" elements."<<endl;

  return vec;
}

std::vector<int> readVectorFromFilec(std::string filename)
{		    
  std::ifstream ifs(filename.c_str(),std::ios::in);	
  
  if (!ifs.is_open())
      throw std::runtime_error("File "+filename+" could not be opened.");

  std::string line;			   
  std::vector<int> vec(0);
  while(!ifs.eof())
    { 
      getline(ifs,line);
       if(line=="")
	break;
      std::stringstream sline(line);
      float val;
      sline>>val;
      vec.push_back(val);
    }
  ifs.close();

  cout<<"Vector "<<filename<<" read, "<<vec.size()<<" elements."<<endl;

  return vec;
}








template <typename T>
void writeVector(std::vector<float> vec,std::string name="")
{
  if(name!="")
    cout<<"Vector "<<name<<" ("<<vec.size()<<"):"<<endl;
  for(unsigned int i=0;i<vec.size();i++)
    cout<<vec[i]<<endl;
  if(name!="")
    cout<<"Vector "<<name<<" end"<<endl;
}

template <typename T>
void writeVector(float* vec,unsigned int length,std::string name="")
{
  if(name!="")
    cout<<"Vector "<<name<<" ("<<length<<"):"<<endl;
  for(unsigned int i=0;i<length;i++)
    cout<<vec[i]<<endl;
  if(name!="")
    cout<<"Vector "<<name<<" end"<<endl;
}


template <typename T>
void writeNonZeroEntries(float* vec, unsigned int length,std::string name="")
{
  if(name!="")
    cout<<"Vector "<<name<<" ("<<length<<"):"<<endl;
  for(unsigned int i=0;i<length;i++)
    if(vec[i])
      cout<<i<<": "<<vec[i]<<endl;
  if(name!="")
    cout<<"Vector "<<name<<" end"<<endl;
}
////  End of IO //////
unsigned int size()
{
	return matrix_size_;
}

void setDoEigenvectors(bool doit) 
{
	do_eigenvectors_=doit;
}
void setDiagonalizationInterval(unsigned int interval) 
{
	diagonalization_interval_=interval;
}
  void setMaxIterations(unsigned int iters) 
  {
	  max_iterations_=iters;
  }
  void setEpsilon(float eps) 
  {
	  epsilon_=eps;
  }
  vector<vector<float> > getEigenvectors()
		  {
	  return eigenvectors_;
		  }
    vector<float> getEigenvalues(){
    	return eigenvalues_;
    }
    void solver(std::vector<float> e, std::vector<float> k)
    {
      n_=e.size()-1;
      
      if(k.size()-1!=n_-1)
        throw std::runtime_error("Jacobi n");

     
      matrix_=new float*[n_];
      for(unsigned int i=0;i<n_;i++)
        {
          matrix_[i]=new float[n_];
          for(unsigned int j=0;j<n_;j++)
    	matrix_[i][j]=0.0;
          matrix_[i][i]=e[i];
          if(i>0)
    	{
    	  matrix_[i][i-1]=k[i-1];
    	  matrix_[i-1][i]=k[i-1];
    	}
        }
     /* for(int i = 0; i < n_; i++) {
       for(int j = 0; j < n_; j++) {
    	cout << matrix_[i][j] << " ";
    	}
    	cout << endl; 
    	fflush(0);
    	}*/
      delete_matrix_=true;
      max_iterations_=100;
    }

    
   //Jacobi
    
    void Doo()
    {
      if(VERBOSITY>1)
      cout<<"Doing Jacobi (size: "<<n_<<")"<<endl;
      float **eigvecs;
      eigvecs=new float *[n_];
      for(unsigned int i=0;i<n_;i++)
        eigvecs[i]=new float[n_];
      int fehler=0;
      //NRjacobi(matrix_,n_,max_iterations_,eigvecs,&fehler);
	  jacobi_method(matrix_,eigvecs,n_);
      if(fehler)
        throw std::runtime_error("jacobi not converged");
      eigenvectors_.clear();
      eigenvalues_.clear();
      for(unsigned int i=0;i<n_;i++)
        {
          std::vector<float> vec(n_);
          for(unsigned int j=0;j<n_;j++)
    	vec[j]=eigvecs[j][i];
          eigenvectors_.push_back(vec);
          eigenvalues_.push_back(matrix_[i][i]);
        }
      for(unsigned int i=0;i<n_;i++)
        delete eigvecs[i];
      delete eigvecs;

    }
    void NRjacobi(float **a, int n, int tmax, float **eigvek, int *fehler)
    {
	
      int  i,j,t,k;
      float  summe,grenze,theta,g,h,tfac,sinus,cosin,speich;
      *fehler=0;
      if(n<2 ) {      
        eigvek[0][0]=1.0;
        return;
      }
      for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
          if(i==j) eigvek[i][i]=1.0; else eigvek[i][j]=0.0; 
        }
      }
      t=0;
      do {
        t++;
        summe=0.0;
        for(i=1;i<n-1;i++) {
          for(j=i+1;j<n;j++) summe+=2*a[i][j]*a[i][j];
        }
      //  cout<<"Helllp"<<endl;
        if(summe==0.0)
        	{
        	return; 
        	}
        else
        {
        	grenze=sqrt(summe)/n/n;
        }
   
        for (i=0;i<n-1;i++) {
          for (j=i+1;j<n;j++) {
            g=100*myabs(a[i][j]);
            if(myabs(a[i][j]) > grenze) {
              h=a[i][i]-a[j][j];
    // This if statement prevents an overflow of theta*theta in the following
    // statement.
              if(myabs(h)+g == myabs(h)) tfac=a[i][j]/h;
              else {
    	    theta=h/2.0/a[i][j];
                tfac=1.0/(myabs(theta)+sqrt(1.0+theta*theta));
                if(theta<0.0)tfac=-tfac;             
              }
              cosin=1.0/sqrt(tfac*tfac+1.0);
              sinus=tfac*cosin;

              for(k=i+1;k<j;k++) {
                speich=a[i][k];
                a[i][k]=cosin*a[i][k]+sinus*a[k][j];
                a[k][j]=cosin*a[k][j]-sinus*speich;
              }
              for(k=j+1;k<n;k++) {
                speich=a[i][k];
                a[i][k]=cosin*a[i][k]+sinus*a[j][k];
                a[j][k]=cosin*a[j][k]-sinus*speich;
              }
              for(k=0;k<i;k++) {
                speich=a[k][i];
                a[k][i]=cosin*a[k][i]+sinus*a[k][j];
                a[k][j]=cosin*a[k][j]-sinus*speich;
              }
              speich=a[i][i];
              a[i][i]=cosin*cosin*a[i][i] + 2*cosin*sinus*a[i][j] +
                      sinus*sinus*a[j][j];
              a[j][j]=cosin*cosin*a[j][j] - 2*cosin*sinus*a[i][j] +
                      sinus*sinus*speich; 
              a[i][j]=0.0;
              for(k=0;k<n;k++) {
                speich=eigvek[k][j];
                eigvek[k][j]=cosin*eigvek[k][j] - sinus*eigvek[k][i];
                eigvek[k][i]=cosin*eigvek[k][i] + sinus*speich;
              }
            }
          }
        }
      } while (t<=tmax);
      *fehler=1;
      cout<<"Jacob its:"<<t<<" tmax"<<tmax<<endl;
      return;

    }


    void sortEigenvalues()
    {
      if(eigenvalues_.size()==0)
        throw std::runtime_error("TridiagonalSolver::sortEigenvalues nothing here.");

    	
      std::vector<unsigned int> positions(0);
    	
      for(unsigned int i=0;i<eigenvalues_.size();i++)
        {
          bool inserted=false;
          for(std::vector<unsigned int>::iterator it=positions.begin();it!=positions.end();it++)
    	{
    	  if(eigenvalues_[*it]>eigenvalues_[i])
    	    {
    	      inserted=true;
    	      positions.insert(it,i);
    	      break;
    	    }
    	}
          if(!inserted)
    	positions.push_back(i);
        }
      std::vector<float> old_vals=eigenvalues_;
      std::vector<std::vector<float> >old_vecs=eigenvectors_;
      for(unsigned int i=0;i<eigenvalues_.size();i++)
        {
          eigenvalues_[i]=old_vals[positions[i]];
          eigenvectors_[i]=old_vecs[positions[i]];
        }

    		
    	
    }
    void calculateEigenvectors(std::vector<float> init_vec,std::vector<std::vector<float> > reduced_eigenvectors)
    {

      float** eigenvectors=new float*[number_eigenvalues_];
      for(unsigned int i=0;i<number_eigenvalues_;i++)
        {
          eigenvectors[i]=new float[size()];
          for(unsigned int j=0;j<size();j++)
    	eigenvectors[i][j]=0.0;
        }

      Initialize(init_vec);
      for(unsigned int step=0;step<e_.size();step++)
        {
          //DoStep();
         
          for(unsigned int i=0;i<number_eigenvalues_;i++)
    	for(unsigned int j=0;j<size();j++)
    	  eigenvectors[i][j]+=act_[j]*reduced_eigenvectors[i][step];
        }
      
      for(unsigned int i=0;i<number_eigenvalues_;i++)
        {
    	  
          std::vector<float> eigenvector(size());
          for(unsigned int j=0;j<size();j++)
    	eigenvector[j]=eigenvectors[i][j];
          
          eigenvectors_.push_back(eigenvector);
    	  
           delete[] eigenvectors[i];
          
        }
      delete[] eigenvectors;
    }

    void jacobi_method ( float ** A, float ** R, int n )
{
// Setting up the eigenvector matrix
for ( int i = 0; i < n; i++ ) {
for ( int j = 0; j < n; j++ ) {
if ( i == j ) {
R[i][j] = 1.0;
} else {
R[i][j] = 0.0;
}
}
}
int k, l;
float epsilon = 1.0e-8;
float max_number_iterations = (float) n * (float) n * (float) n;
int iterations = 0;
float max_offdiag = maxoffdiag ( A, &k, &l, n );
while ( fabs(max_offdiag) > epsilon && (float) iterations < max_number_iterations ) {
max_offdiag = maxoffdiag ( A, &k, &l, n );
rotate ( A, R, k, l, n );
iterations++;
}
std::cout << "Number of iterations: " << iterations << "\n";
return;
}

float maxoffdiag ( float ** A, int * k, int * l, int n )
{
float max = 0.0;
for ( int i = 0; i < n; i++ ) {
for ( int j = i + 1; j < n; j++ ) {
if ( fabs(A[i][j]) > max ) {
max = fabs(A[i][j]);
*l = i;
*k = j;
}
}
}
return max;
}
// Function to find the values of cos and sin
void rotate ( float ** A, float ** R, int k, int l, int n )
{
float s, c;
if ( A[k][l] != 0.0 ) {
float t, tau;
tau = (A[l][l] - A[k][k])/(2*A[k][l]);
if ( tau > 0 ) {
t = 1.0/(tau + sqrt(1.0 + tau*tau));
} else {
t = -1.0/( -tau + sqrt(1.0 + tau*tau));
}
c = 1/sqrt(1+t*t);
s=c*t;
} else {
c = 1.0;
s = 0.0;
}
float a_kk, a_ll, a_ik, a_il, r_ik, r_il;
a_kk = A[k][k];
a_ll = A[l][l];
// changing the matrix elements with indices k and l
A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll;
A[k][l] = 0.0; // hard-coding of the zeros
A[l][k] = 0.0;
// and then we change the remaining elements
for ( int i = 0; i < n; i++ ) {
if ( i != k && i != l ) {
a_ik = A[i][k];
a_il = A[i][l];
A[i][k] = c*a_ik - s*a_il;
A[k][i] = A[i][k];
A[i][l] = c*a_il + s*a_ik;
A[l][i] = A[i][l];
}
// Finally, we compute the new eigenvectors
r_ik = R[i][k];
r_il = R[i][l];
R[i][k] = c*r_ik - s*r_il;
R[i][l] = c*r_il + s*r_ik;
}
return;
}
    
    
    
    
    
    
