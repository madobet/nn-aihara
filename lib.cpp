#include <iostream>
#include "lib.h"

namespace ntoolkit{

  std::string& IO::ReplaceAllStr(std::string& str, const std::string& old_str, const std::string& new_str){
      while(true) {
          std::string::size_type pos(0);
          if( (pos = str.find(old_str)) != std::string::npos ) str.replace(pos,old_str.length(),new_str);
          else break;
      }
      return str;
  }

  template <typename T>
  std::vector< std::vector<T> > IO::ReadMatFromFile(const std::string &filepath){
      std::ifstream wfs(filepath);
      if(!wfs.is_open()) {
          wfs.close();
          throw "weight file not found" ;
      }
      int net_size = GetMatDimFromFile(filepath);
      if(net_size < 0){
          throw "can not analysis matrix dimention from file";
      }

      // Start reading from file
      std::vector< std::vector<T> > res; res.clear();
      for (int i = 0; i < net_size; i++)
      {
          std::vector<T> tmp_vec;
          tmp_vec.clear();
          for (int j = 0; j < net_size; j++)
          {
              double temp;
              wfs >> temp;
              tmp_vec.push_back(temp);
          }
          res.push_back(tmp_vec);
      }
      wfs.close();
      return res;
  }

  template <typename T>
  std::vector< std::vector<T> > IO::ReadMatFromFile(const std::string &filepath, T threshold){
      std::ifstream wfs(filepath);
      if(!wfs.is_open()) {
          wfs.close();
          throw "weight file not found" ;
      }
      int net_size = GetMatDimFromFile(filepath);
      if(net_size < 0){
          throw "can not analysis matrix dimention from file";
      }

      // Start reading from file
      std::vector< std::vector<T> > res; res.clear();
      for (int i = 0; i < net_size; i++)
      {
          std::vector<T> tmp_vec;
          tmp_vec.clear();
          for (int j = 0; j < net_size; j++)
          {
              double temp;
              wfs >> temp;
              if(abs(temp) > threshold)   tmp_vec.push_back(temp);
              else                        tmp_vec.push_back((double)0.0);
          }
          res.push_back(tmp_vec);
      }
      wfs.close();
      return res;
  }

  template <typename T>
  bool IO::WriteMat2File(std::vector<std::vector<T> > &m, std::string &fpath){
      std::ofstream ofs(fpath);
      for(const std::vector<T> &b : m ){
          for(const T &a : b){
              ofs << a << " ";
          }
          ofs << std::endl;
      }
      ofs.close();
  }

  int IO::GetMatDimFromFile(const std::string &filepath){
      std::ifstream infile(filepath);
      std::string s;
      int n, res0, res1;
      res0 = 0;
      getline(infile, s);
      std::istringstream iss(s);
      while(iss >> n)
          res0++;
      while(getline(infile, s)){
          res1 = 0;
          iss.str(s);
          while(iss >> n)
              res1++;
          if(res1!=res0){
              return -1;
          }
          res0 = res1;
      }
      infile.close();
      return res1;
  }

  template <typename T>
  bool IO::SaveVec2File(std::vector<T> &array, std::string &filename){
      std::ofstream ofs(filename);
      for(const T &b : array ){
          ofs << b << " ";
      }
      ofs << std::endl;
      ofs.close();
      return true;
  }

  template <typename T>
  std::vector<T> IO::ExtOneColFromMat(std::vector< std::vector<T> > &m, int cid){
      std::vector<T> res;
      if(m.empty())
          return std::vector<T>();
      int time_step = m.size();
      for(int i=0; i < time_step; i++){
          res.push_back(m.at(i).at(cid));
      }
      return res;
  }

  const double Math::pi=3.14159265358979323846;
  const double Math::eps=1e-6;
  int Math::RandNum(std::pair<int, int> &lim_lo2up){
      static std::default_random_engine e;
      e.seed(time(0)); //!< @note 静态变量和设置time函数用作种子，双保险，保证随机性
      static std::uniform_int_distribution<int> u(lim_lo2up.first,lim_lo2up.second);
      return u(e);
  }
  double Math::RandNum(std::pair<double, double> &lim_lo2up){
      static std::default_random_engine e;
      e.seed(time(0)); //!< @note 静态变量和设置time函数用作种子，双保险，保证随机性
      static std::uniform_real_distribution<double> u(lim_lo2up.first,lim_lo2up.second);
      return u(e);
  }

  template <typename T>
  std::vector<T> Math::RandVec(T size, std::pair<T, T> &lim_lo2up){
      std::vector<T> res; res.clear();
      for(int i = 0; i < size; i++) res.push_back(RandNum(lim_lo2up));
      return res;
  }

  template <typename T>
  std::vector< std::vector<T> > Math::RandMat(std::pair<T, T> &m_mul_n, std::pair<T, T> &lim_lo2up){
      std::vector<T> res; res.clear();
      for(int i = 0; i < m_mul_n.first; i++){
          std::vector<T> row_vec; row_vec.clear();
          for(int j = 0; j < m_mul_n.second; j++) row_vec.push_back(RandNum(lim_lo2up));
          res.push_back(row_vec);
      }
      return res;
  }

  double Math::CubicEquRealRT(double a, double b, double c, double d){
      double delta = pow(36*a*b*c - 8*pow(b,3) - 108*(pow(a,2)*d),2) \
         + pow(12*a*c - 4*pow(b,2) ,3);
      if(delta > 0);
      double root = 0;  //!< @todo unfinished
      return root;
  }

  double Math::RKMethod() {
      /// @todo 4th R-K Method
      return 0;
  }

  double Math::NewtonCubicEquRT(double x0, double err, double a, double b, double c, double d){
      double xn, xn1, dy, y;
      if( (3*a*pow(x0,2) + 2*b*x0 + c) == 0 ) xn = x0 + 1;
      else xn = x0;

      dy = 3*a*pow(xn,2) + 2*b*xn + c;
      y = a*pow(xn,3) + b*pow(xn,2) + c*xn + d;
      xn1 = xn - y/dy;
      while(fabs(xn1-xn) > err){
          xn = xn1;
          dy = 3*a*pow(xn,2) + 2*b*xn + c;
          y = a*pow(xn,3) + b*pow(xn,2) + c*xn + d;
          xn1 = xn - y/dy;
      }
      return xn1;
  }

  double Math::CalExp(char *expression )
  {
      char temp[max_buffer_size]={0};
      Trans2Inverse(expression, temp);
      return CompValue(temp);
  }

  double Math::CalFunc(std::string &func, const std::string &arg, const std::string &arg_value) {
      IO::ReplaceAllStr(func,"sin","s");
      IO::ReplaceAllStr(func,"cos","c");
      IO::ReplaceAllStr(func,"tan","t");
      IO::ReplaceAllStr(func,"sqrt","v");
      IO::ReplaceAllStr(func,"pow","^");
      IO::ReplaceAllStr(func,"log","L");
      IO::ReplaceAllStr(func,"ln","l");
      IO::ReplaceAllStr(func,"exp","e");
      // Replace all long string for CalExp

      IO::ReplaceAllStr(func,arg,arg_value);
      char * ex = new char[func.size()+1]; // end with 0 so +1
      ex[func.size()]=0;
      double result = CalExp(strcpy(ex,func.c_str()));
      delete []ex;

      return result;
  }


  void Math::Trans2Inverse(char *a, char *b)
  {
      char stock[max_buffer_size]={0};

      int top=0;
      int len=0;
      int i=0;
      int j=0;

      top = -1;
      j = -1;
      len = strlen(a);

      for ( i=0; i<len; i++ )
      {
          switch( a[i] )
          {
              case '(':
                  stock[++top] = '(';
                  break;

              case '+':
              case '-':
                  while( top>=0 && stock[top]!='(' )
                  {
                      b[++j] = stock[top--];
                  }
                  stock[++top] = ' ';
                  stock[++top] = a[i];
                  break;

              case '*':
              case '/':
                  while( top>=0 && stock[top]!='(' && stock[top]!='+' && stock[top]!='-' )
                  {
                      b[++j] = stock[top--];
                  }
                  stock[++top] = ' ';
                  stock[++top] = a[i];
                  break;

              case 's':
              case 'c':
              case 't':
                  while( top>=0 && stock[top]!='(' && stock[top]!='+' && stock[top]!='-' && stock[top]!='*' && stock[top]!='/' )
                  {
                      b[++j] = stock[top--];
                  }
                  stock[++top] = ' ';
                  stock[++top] = a[i];
                  break;

              case 'v':
              case '^':
                  while( top>=0 && stock[top]!='(' && stock[top]!='+' && stock[top]!='-' && stock[top]!='*' && stock[top]!='/' && stock[top]!='s' && stock[top]!='c' && stock[top]!='t' )
                  {
                      b[++j] = stock[top--];
                  }
                  stock[++top] = ' ';
                  stock[++top] = a[i];
                  break;

              case 'L':
              case 'l':
              case 'e':
                  while( top>=0 && stock[top]!='(' && stock[top]!='+' && stock[top]!='-' && stock[top]!='*' && stock[top]!='/' && stock[top]!='s' && stock[top]!='c' && stock[top]!='t' && stock[top]!='v' && stock[top]!='^' )
                  {
                      b[++j] = stock[top--];
                  }
                  stock[++top] = ' ';
                  stock[++top] = a[i];
                  break;

              case')':
                  while( stock[top]!='(' )
                  {
                      b[++j] = stock[top--];
                  }
                  top--;
                  break;

              default:
                  b[++j] = a[i];
                  if( i == len-1 || a[i+1]<'0' || a[i+1]>'9' )
                  {
                      if ( a[i+1] != '.' )
                      {
                          b[++j] = ' ';
                      }
                  }
                  break;
          }
      }

      while ( top>=0 )
      {
          b[++j] = stock[top--];
      }

      b[++j] = '\0';
  }
  
  double Math::CompValue(char *expression)
  {
      int top=0;
      int len=0;
      int i=0;
      int c=0;

      double sum=0;
      double digit[max_buffer_size]={0};

      char str_num_temp[max_buffer_size]={0};

      top = -1;
      len = strlen(expression);

      for ( i=0; i<len; i++ )
      {
          switch( expression[i] )
          {
              case ' ':
                  break;

              case '+':
                  sum = digit[top] + digit[top-1];
                  digit[--top] = sum;
                  break;

              case '-':
                  sum = digit[top-1] - digit[top];
                  digit[--top] = sum;
                  break;

              case '*':
                  sum = digit[top] * digit[top-1];
                  digit[--top] = sum;
                  break;

              case '/':
                  sum = digit[top-1] / digit[top];
                  digit[--top] = sum;
                  break;

              case 's':
                  sum = sin(Deg2Rad(digit[top]) );
                  digit[top] = sum;
                  break;

              case 'c':
                  sum = cos(Deg2Rad(digit[top]) );
                  digit[top] = sum;
                  break;

              case 't':
                  sum = tan(Deg2Rad(digit[top]) );
                  digit[top] = sum;
                  break;

              case 'v':
                  sum = sqrt( digit[top] );
                  digit[top] = sum;
                  break;

              case '^':
                  sum = pow( digit[top-1], digit[top] );
                  digit[--top] = sum;
                  break;

              case 'L':
                  sum = log10( digit[top] );
                  digit[top] = sum;
                  break;

              case 'l':
                  sum = log( digit[top] );
                  digit[top] = sum;
                  break;

              case 'e':
                  sum = exp( digit[top] );
                  digit[top] = sum;
                  break;

              default:
                  c = 0;
                  memset( str_num_temp, 0, sizeof(str_num_temp) );
                  while( expression[i]>='0' && expression[i]<='9' || expression[i] == '.' )
                  {
                      str_num_temp[c++] = expression[i];
                      i++;
                  }
                  str_num_temp[c] = '\0';
                  digit[++top] = strtod(str_num_temp, nullptr);
                  break;
          }
      }

      return digit[0];
  }

  Graph::Graph(){
      igraph_empty(&graf, 0, IGRAPH_UNDIRECTED);
      igraph_matrix_init(&wmat, 0, 0);
      igraph_matrix_null(&wmat);
  }

  Graph::Graph(const std::vector<std::vector<double>> &m, bool directed){
      int s;
      if(m.empty()){
          std::cerr << "Error:Matrix empty , set graph to 0" << std::endl;
          s = 0;
      }else{
          s = (int)m.size();
      }
      if(directed){
          igraph_empty(&graf, s, IGRAPH_DIRECTED);
      }else{
          igraph_empty(&graf, s, IGRAPH_UNDIRECTED);
      }
      igraph_matrix_init(&wmat, s, s);
      for(int i=0; i<s; i++) //!< @remark i表示行标，j表示列标
          for(int j=0; j<s; j++){
              igraph_matrix_set(&wmat, i, j, m.at(i).at(j));
          }//!< @note 如果igraph_matrix_set线程安全，那么可以并行
  }

  Graph::Graph(std::vector<int> &scale_vec, std::vector<double> &weight) {
      igraph_vector_t dim_vec;
      Vec2igVec(scale_vec, dim_vec);
      igraph_vector_destroy(&dim_vec);
  }

  Graph::Graph(const Graph & G)
  {
      /*
       * 此函数深度复制图形对象以创建它的精确副本。
       * 新的副本应该igraph_destroy()在不再需要时通过调用来销毁 。
       * 您还可以通过简单地使用标准的赋值运算符创建一个图形的浅拷贝，
       * 但要小心，不要不 破坏浅副本。
       * 为避免此错误，建议不要创建浅拷贝。
       * */
      igraph_copy(&graf, &G.graf);

      igraph_matrix_copy(&wmat, &G.wmat);
  }

  Graph& Graph::operator =(const Graph &G){
      if(this==&G)
          return *this;
      igraph_destroy(&(this->graf));
      igraph_copy(&(this->graf), &G.graf);
      igraph_matrix_destroy(&(this->wmat));
      igraph_matrix_copy(&wmat, &G.wmat);
      return *this;
  }

  Graph::~Graph() {
      igraph_matrix_destroy(&wmat);
      /*
       * 应该为每个图形对象调用此函数一次。
       * 这个函数使所有迭代器失效（当然），但是图形的迭代器应该在图形本身之前销毁。
       */
      igraph_destroy(&graf);
  }

  bool Graph::ConvertDirectedness(bool directed) {
      if(directed){
          igraph_to_directed(&graf, IGRAPH_TO_DIRECTED_MUTUAL);
      }else{
          igraph_to_undirected(&graf, IGRAPH_TO_UNDIRECTED_MUTUAL, nullptr);
          //!< 会丢失原图属性？
      }
      return this->IsDirected();
  }

  std::pair<int,int> Graph::EdgeEndpoints(int edge_id) const {
      igraph_integer_t s_id;
      igraph_integer_t e_id;
      igraph_edge(&graf, edge_id, &s_id, &e_id);
      return std::make_pair((int)s_id,(int)e_id);
  }

  int Graph::EdgeBetween(const std::pair<int, int> &vid_pair ) const {
      igraph_integer_t eid;
      if(igraph_get_eid(&graf, &eid, vid_pair.first, vid_pair.second, igraph_is_directed(&graf), true)==IGRAPH_SUCCESS)
          return eid;
      return invalid;
  }

  std::vector<int>
  Graph::EdgesAlongPath(const std::vector<std::pair<int, int>> &vid_pairs, const std::vector<int> &vid_along_path) const {
      igraph_vector_t pairs;
      igraph_vector_t path;
      igraph_vector_t e_ids;
      igraph_vector_t *pt_pairs=&pairs;
      igraph_vector_t *pt_path=&path;
      std::vector<int> res;
      res.clear();
      if(vid_pairs.empty()){
          pt_pairs = nullptr;
      }else{
          FlatVecPairs2igVec(vid_pairs, pairs);
      }
      if(vid_along_path.empty()){
          pt_path = nullptr;
      }else{
          Vec2igVec(vid_along_path, path);
      }
      if(igraph_get_eids_multi(&graf, &e_ids, pt_pairs, pt_path, igraph_is_directed(&graf), true)==IGRAPH_SUCCESS){
          res= IgVec2intVec<int>(e_ids);
      }else{
          res.clear();
      }
      igraph_vector_destroy(&e_ids);
      igraph_vector_destroy(&path);
      igraph_vector_destroy(&pairs);
      return res;
  }

  std::vector<int> Graph::Degree() const {
      // if(igraph_matrix_empty(&wmat))???没看懂这个函数，不敢用
      std::vector<int> degree;
      igraph_vector_t res;
      igraph_vector_init(&res,igraph_vcount(&graf));
      igraph_degree(&graf, &res, igraph_vss_all(), IGRAPH_ALL, false);
      degree = IgVec2intVec<int>(res);
      igraph_vector_destroy(&res);
      return degree;
  }

  std::vector<int> Graph::AloneVertices() const {
      std::vector<int> vid;
      vid.clear();
      std::vector<int> all = Degree();
      for(unsigned int i=0; i<igraph_vcount(&graf); i++){
          if(all.at(i)==0)
              vid.push_back(i);
      }
      return vid;
  }

  bool Graph::IsStrongConnected() {
      igraph_bool_t res;
      igraph_is_connected(&graf, &res, IGRAPH_STRONG);
      return res?true:false;
  }

  bool Graph::IsWeakConnected() {
      igraph_bool_t res;
      igraph_is_connected(&graf, &res, IGRAPH_WEAK);
      return res?true:false;
  }

  bool Graph::AddEdges(std::vector<std::pair<std::pair<int, int> , double>>  &list) {
      igraph_vector_t vec;
      std::vector<std::pair<int, int>> vlist = ExtVecPairs1stVec(list);
      FlatVecPairs2igVec(vlist, vec);
      switch(igraph_add_edges(&graf, &vec, nullptr)){
          case IGRAPH_EINVEVECTOR:
              std::cerr << "Error: invalid (odd) edges vector length" << std::endl;
              igraph_vector_destroy(&vec);
              return false;
          case IGRAPH_EINVVID:
              std::cerr << "Error: invalid vertex id in edges vector" << std::endl;
              igraph_vector_destroy(&vec);
              return false;
          default:
              igraph_vector_destroy(&vec);
              break;
      }
      for( const auto &it : list ){
          UpdateWMat(it.first.first, it.first.second, it.second);
      }
      return true;
  }

  bool Graph::AddVertices(int nv) {
      if(igraph_add_vertices(&graf, nv, nullptr) == IGRAPH_EINVAL){
          std::cerr << "Error invalid number of new vertices" <<std::endl;
          return false;
      }

      return true;
  }

  bool Graph::RandAddEdges(std::vector<double> &wlist) {
      igraph_vector_t evid;
      igraph_rng_seed(igraph_rng_default(), 9323);
      igraph_vector_init(&evid, wlist.size()*2);
      for (int i=0; i<wlist.size(); i++) {
          std::default_random_engine e;
          e.seed(time(nullptr));
          std::uniform_int_distribution<int> u(0,(int)igraph_vcount(&graf));
          int fst = u(e);
          int snd = u(e);
          while(fst == snd){
              snd = u(e);
          } //!< 防止产生自环
          VECTOR(evid)[2*i] = fst;
          VECTOR(evid)[2*i+1] = snd;
      }
      switch(igraph_add_edges(&graf, &evid, nullptr)){
          case IGRAPH_EINVEVECTOR:
              std::cerr << "Error: invalid (odd) edges vector length" << std::endl;
              igraph_vector_destroy(&evid);
              return false;
          case IGRAPH_EINVVID:
              std::cerr << "Error: invalid vertex id in edges vector" << std::endl;
              igraph_vector_destroy(&evid);
              return false;
          default:
              break;
      }
      for (int i=0; i<wlist.size(); i++){
          UpdateWMat(VECTOR(evid)[2*i], VECTOR(evid)[2*i+1], wlist.at(i));
      }
      igraph_vector_destroy(&evid);
      return true;
  }

  bool Graph::DelEdges(std::vector<std::pair<int, int>> &vlist) {
      igraph_vector_t v;
      FlatVecPairs2igVec(vlist, v);
      igraph_es_t egs;
      igraph_es_pairs(&egs, &v, igraph_is_directed(&graf));
      igraph_delete_edges(&graf, egs);
      igraph_es_destroy(&egs);
      igraph_vector_destroy(&v);
      for( const auto &it : vlist ){
          UpdateWMat(it.first, it.second, 0);
      }
      return true;
  }

  bool Graph::DelVertices(std::vector<int> &vid) {
      igraph_vector_t v;
      Vec2igVec(vid, v);
      if(igraph_delete_vertices(&graf, igraph_vss_vector(&v)) == IGRAPH_EINVVID){
          std::cerr << "Error: invalid vertex id" <<std::endl;
          igraph_vector_destroy(&v);
          return false;
      }
      return true;
  }

  double Graph::UnweightShortPaths(const std::pair<int, int> &vid_pair) const {
      igraph_matrix_t mat;
      igraph_vector_t from;
      igraph_vector_t to;
      igraph_matrix_init(&mat, 1, 1);
      igraph_vector_init(&from, 1);
      igraph_vector_init(&to, 1);
      VECTOR(from)[0] = vid_pair.first;
      VECTOR(to)[0] = vid_pair.second;
      double res = invalid;
      switch(igraph_shortest_paths(&graf, &mat, igraph_vss_vector(&from), igraph_vss_vector(&to), IGRAPH_OUT)){
          case IGRAPH_ENOMEM:
              std::cerr << "ERROR: not enough memory for temporary data" << std::endl;
              break;
          case IGRAPH_EINVVID:
              std::cerr << "ERROR: invalid vertex id passed" << std::endl;
              break;
          case IGRAPH_EINVMODE:
              std::cerr << "ERROR: invalid mode argument" << std::endl;
              break;
          default:
              res = MATRIX(mat, 1, 1);
              break;
      }
      igraph_vector_destroy(&to);
      igraph_vector_destroy(&from);
      igraph_matrix_destroy(&mat);
      return res;
  }

  double Graph::AverPathLength() const {
      double res;
      igraph_average_path_length(&graf, &res, true, false);
      return res;
  }

  std::vector<double> Graph::PathLengthDist() const {
      std::vector<double> res;
      res.clear();
      double unreach; // 无法相互连通的对的数量
      igraph_vector_t vec;
      igraph_vector_init(&vec, 1);
      igraph_path_length_hist(&graf, &vec, &unreach, true);
      res.push_back(unreach); // 把无法连通，也即长度为零的放在零位置上
      for(int i=1; i<igraph_vector_size(&vec); i++){
          res.push_back(VECTOR(vec)[i-1]);
      }
      res = IgVec2intVec<double>(vec);
      igraph_vector_destroy(&vec);
      return res;
  }

  double Graph::Diameter() const {
      return 0;
  }

  std::vector<std::vector<double>>
  Graph::UnweightShortPaths(const std::pair<std::vector<int>, std::vector<int>> &pair_vids) const {
      igraph_matrix_t mat;
      igraph_vector_t from;
      igraph_vector_t to;
      igraph_matrix_init(&mat, pair_vids.first.size(), pair_vids.second.size());
      Vec2igVec(pair_vids.first,from);
      Vec2igVec(pair_vids.second,to);
      std::vector<std::vector<double>> res;
      res.clear();
      switch(igraph_shortest_paths(&graf, &mat, igraph_vss_vector(&from), igraph_vss_vector(&to), IGRAPH_OUT)){
          case IGRAPH_ENOMEM:
              std::cerr << "ERROR: not enough memory for temporary data" << std::endl;
              break;
          case IGRAPH_EINVVID:
              std::cerr << "ERROR: invalid vertex id passed" << std::endl;
              break;
          case IGRAPH_EINVMODE:
              std::cerr << "ERROR: invalid mode argument" << std::endl;
              break;
          default:
              res = IgMatrix2VecVecs<double>(mat);
              break;
      }
      igraph_vector_destroy(&to);
      igraph_vector_destroy(&from);
      igraph_matrix_destroy(&mat);
      return res;
  }

  std::vector<std::vector<double>>
  Graph::UnweightShortPaths() const {
      igraph_matrix_t mat;
      igraph_matrix_init(&mat, igraph_vcount(&graf), igraph_vcount(&graf));
      std::vector<std::vector<double>> res;
      res.clear();
      switch(igraph_shortest_paths(&graf, &mat, igraph_vss_all(), igraph_vss_all(), IGRAPH_OUT)){
          case IGRAPH_ENOMEM:
              std::cerr << "ERROR: not enough memory for temporary data" << std::endl;
              break;
          case IGRAPH_EINVVID:
              std::cerr << "ERROR: invalid vertex id passed" << std::endl;
              break;
          case IGRAPH_EINVMODE:
              std::cerr << "ERROR: invalid mode argument" << std::endl;
              break;
          default:
              res = IgMatrix2VecVecs<double>(mat);
              break;
      }
      igraph_matrix_destroy(&mat);
      return res;
  }

  std::vector<int> Graph::ArticulationPoints() const {
      std::vector<int> res;
      res.clear();
      igraph_vector_t vec;
      igraph_vector_init(&vec, 0);
      igraph_articulation_points(&graf, &vec);
      res = IgVec2intVec<int>(vec);
      igraph_vector_destroy(&vec);
      return res;
  }

  int Graph::MaxDegree() const {
      int res = 0;
      switch(igraph_maxdegree(&graf, &res, igraph_vss_all(), IGRAPH_ALL, false)){
          case IGRAPH_EINVVID:
              std::cerr << "invalid vertex id" << std::endl;
              break;
          case IGRAPH_EINVMODE:
              std::cerr << "invalid mode argument" << std::endl;
              break;
          default:
              break;
      }
      return res;
  }

  std::vector<int> Graph::MaxDegreePoints() const {
      std::vector<int> res;
      res.clear();
      std::vector<int> de = Degree();
      int max_de = MaxDegree();
      for(int i=0; i<igraph_vcount(&graf); i++){
          if(de.at(i)==max_de){
              res.push_back(i);
          }
      }
      return res;
  }

  double Graph::ClusterCoeff() const {
      double res;
      switch(igraph_transitivity_undirected(&graf,&res,IGRAPH_TRANSITIVITY_ZERO)){
          case IGRAPH_ENOMEM:
              std::cerr << "ERROR: not enough memory for temporary data" << std::endl;
              break;
          default:
              break;
      };
      return res;
  }

  std::vector<std::vector<int>> Graph::MiniumSizeSeps() const {
      std::vector<std::vector<int>> res;
      res.clear();
      igraph_vector_ptr_t seps;
      igraph_vector_ptr_init(&seps, 0);
      igraph_minimum_size_separators(&graf,&seps);
      for(int i=0; i<igraph_vector_ptr_size(&seps); i++){
          igraph_vector_t * v = (igraph_vector_t*)VECTOR(seps)[i];
          res.push_back(IgVec2intVec<int>(*v));
          igraph_vector_destroy(v);
          igraph_free(v);
      }
      igraph_vector_ptr_destroy(&seps);
      return res;
  }

  template <typename T>
  void Graph::Vec2igVec(const std::vector<T> &std_vec, igraph_vector_t &igraph_vec) const {
      igraph_vector_init(&igraph_vec, std_vec.size());
      if(std_vec.empty()){
          igraph_vector_destroy(&igraph_vec);
          return;
      }
      for(int i=0; i<igraph_vector_size(&igraph_vec); i++){
          VECTOR(igraph_vec)[i] = std_vec.at(i);
      }
  }

  template <typename T>
  std::vector<T> Graph::IgVec2intVec(const igraph_vector_t &igraph_vec) const {
      std::vector<T> res;
      res.clear();
      for (int i=0; i<igraph_vector_size(&igraph_vec); i++) {
          T element = VECTOR(igraph_vec)[i];
          res.push_back(element);
      }
      return res;
  }

  template <typename T>
  void Graph::FlatVecPairs2igVec(const std::vector<std::pair<T, T>> &src, igraph_vector_t &dst) const {
      igraph_vector_init(&dst, src.size()*2);
      for(int i=0; i<src.size(); i++){
          VECTOR(dst)[2*i] = src.at(i).first;
          VECTOR(dst)[2*i+1] = src.at(i).second;
      }
  }
  template <typename T>
  void Graph::DivVecPairs2igVecs(const std::vector<std::pair<T, T>> &src, igraph_vector_t &dst1,
                                 igraph_vector_t &dst2) const {
      igraph_vector_init(&dst1, src.size());
      igraph_vector_init(&dst2, src.size());
      for(int i=0; i<src.size(); i++){
          VECTOR(dst1)[i] = src.at(i).first;
          VECTOR(dst2)[i] = src.at(i).second;
      }
  }

  template<typename T, typename U>
  void Graph::DivPairVecs2igVecs(const std::pair<std::vector<T>, std::vector<U>> &src, igraph_vector_t &dst1,
                                 igraph_vector_t &dst2) const {
      igraph_vector_init(&dst1, src.first.size());
      for(int i=0; i<src.first.size(); i++){
          VECTOR(dst1)[i] = src.first.at(i);
      }
      igraph_vector_init(&dst2, src.second.size());
      for(int i=0; i<src.second.size(); i++){
          VECTOR(dst2)[i] = src.second.at(i);
      }
  }

  template<typename T>
  std::vector<std::vector<T>> Graph::IgMatrix2VecVecs(const igraph_matrix_t &src) const {
      std::vector<std::vector<T>> res;
      res.clear();
      int nrow = igraph_matrix_nrow(&src);
      int ncol = igraph_matrix_ncol(&src);
      T* to;
      to = new T[nrow*ncol];
      igraph_matrix_copy_to(&src, to); // 是按列复制进去的
      for(int i=0; i<nrow; i++){
          std::vector<T> row;
          row.clear();
          for(int j=0; j<ncol; j++){
              row.push_back(to[i+j*nrow]);
          }
          res.push_back(row);
      }
      delete[] to;
      return res;
  }

  template<typename T, typename U>
  std::vector<T> Graph::ExtVecPairs1stVec(const std::vector<std::pair<T, U>> &src) const {
      std::vector<T> res;
      res.clear();
      for( const std::pair<T,U> &it : src ){
          res.push_back(it.first);
      }
      return res;
  }

  template<typename T, typename U>
  std::vector<U> Graph::ExtVecPairs2ndVec(const std::vector<std::pair<T, U>> &src) const {
      std::vector<U> res;
      res.clear();
      for( const std::pair<T,U> &it : src ){
          res.push_back(it.first);
      }
      return res;
  }
}   // namespace ntoolkit
