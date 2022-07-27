#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

List MaxStats(NumericMatrix X,String optimization, String method, bool circular, int min_seg);
List MaxStatsL1(NumericMatrix X, String method, int min_seg);
List MaxStatsL2(NumericMatrix X, String method, int min_seg);
List BinarySegmentationL1(NumericMatrix X, int min_seg);
List BinarySegmentationL2(NumericMatrix X, int min_seg);
double EvalStatL2(NumericMatrix S, int s, int t);
List OneStepSearchL2(NumericMatrix S, int s, String method, int min_seg);
IntegerVector StoppingBoundary(int nr_perms, double alpha, double eta_star);
bool PermTest(NumericMatrix X, IntegerVector boundary, double candidate_stat,String optimization, String method, double alpha, int nr_perm, bool circular, int min_seg);
double pPerm(NumericMatrix X, IntegerVector boundary, double cand_stat, String optimization, String method, double alpha,
              int nr_perms, bool circular, int min_seg);
IntegerVector ChangePoints(NumericMatrix X, IntegerVector boundary, String optimization, String method, double alpha, int nr_perms, double eta, int min_seg);


IntegerVector vecpow(const int base, IntegerVector exp) {
  // [base^exp_1,..., base^exp_n]
  for(int i=0; i<exp.length(); i++)
    exp[i] = pow(base, exp(i));
  return exp;
}

IntegerVector seqnn(int a, int b){
  // integer sequence non-negative, i.e. empty if b<a
  IntegerVector res;
  if(a<=b)
    res = seq(a,b);
  return res;
}

void merge_vecs(IntegerVector& x, IntegerVector y){
  // merge x and y into x and sort
  int n = y.length();
  for(int i=0; i<n; i++){
    x.push_back(y(i));
  }
  sort(x.begin(),x.end());
}

class BinTree{
// Tree structure used for L1 loss
private:
  NumericMatrix X;
  int T;
  int n;
  int k;
  IntegerVector boundaries;
  IntegerMatrix tree;
  int len_prev = 0;
  IntegerMatrix tree_prev;
  NumericVector loss;
  NumericVector loss_prev;
  IntegerMatrix tree_to_orig;
  IntegerMatrix orig_to_tree;

public:
  BinTree(NumericMatrix X_in){
    // constructor
    X = X_in;
    T = X.ncol();
    n = X.nrow();
    k = ceil(log2(T));

    // tree[boundaries[i:i+1]] contains information of tree height i
    boundaries = vecpow(2,k-seq(0,k-1));
    partial_sum(boundaries.begin(), boundaries.end(), boundaries.begin());
    boundaries.push_front(0);

    // Tree structure to save and evaluate efficiently which nodes are active
    IntegerMatrix _tree(n,2*pow(2,k)-1);
    tree = _tree;
    if(pow(2,k)>T){
      for(int t:seqnn(T, pow(2,k)-1))
        tree(_,t) = rep(R_NaN,n);
    }
    tree_prev = tree;

    // tree_to_orig[0] is the location of the smallest element in X
    IntegerMatrix _tree_to_orig(n,T);
    tree_to_orig = _tree_to_orig;
    for(int i:seqnn(0,n-1)){
      IntegerVector tree_to_orig_i = seq(0,T-1);
      sort(tree_to_orig_i.begin(), tree_to_orig_i.end(),
           [&](const int& a, const int& b) {
             return (X(i,a) < X(i,b));
           }
      );
      tree_to_orig(i,_) = tree_to_orig_i;
    }

    // orig_to_tree[0] is the location of the first element of X in the ordered tree
    IntegerMatrix _orig_to_tree(n,T);
    orig_to_tree = _orig_to_tree;
    for(int i:seqnn(0,n-1)){
      IntegerVector orig_to_tree_i = seq(0,T-1);
      sort(orig_to_tree_i.begin(), orig_to_tree_i.end(),
           [&](const int& a, const int& b) {
             return (tree_to_orig(i,a) < tree_to_orig(i,b));
           }
      );
      orig_to_tree(i,_) = orig_to_tree_i;
    }

    // loss = [loss of inactive nodes, loss of active nodes, total loss]
    double loss0 = 0;
    for(int i:seqnn(0,n-1)){
      double med0;
      if(T % 2 == 0)
        med0 = (X(i, tree_to_orig(i, T/2-1))+ X(i, tree_to_orig(i, T/2)))/2;
      else
        med0 = X(i, tree_to_orig(i, (T-1)/2));
      loss0 += sum(abs(X(i,_)-med0));
    }

    loss = NumericVector::create(loss0, 0 , loss0);
    loss_prev = loss;
  }

public:
  void Flip(IntegerVector pts);
  // void Flip2(int pt1, int pt2);
  double EvalStatL1();
  List MaxStatsL1(String method, int min_seg);
  List BinarySegmentationL1(int min_seg);

private:
  void ResetTree(bool prev);
  int NodesOf(bool type);
  NumericVector FindMedian(bool type);
  IntegerVector mthElement(bool type, int m);
  List OneStepSearchL1(int len);
};

void BinTree::ResetTree(bool prev){
  // Sets the tree either to empty (prev=0) or to the start of the previous
  // diagonal (prev=1)
  if(prev==0){
    IntegerMatrix _tree(n,2*pow(2,k)-1);
    tree = _tree;
    if(pow(2,k)>T){
      for(int t:seqnn(T, pow(2,k)-1))
        tree(_,t) = rep(R_NaN,n);
    }
    loss = NumericVector::create(loss(2), 0, loss(2));
  }
  else{
    tree = tree_prev;
    loss = loss_prev;
  }

}

int BinTree::NodesOf(bool type){
  // Since the number of nodes is the same in all dimensions at all times,
  // we return the number of active nodes in the first dimension
  // For type == 0, the number of inactive nodes are returned
  if(type==1)
    return(tree(0, boundaries(k)));
  else
    return(T-tree(0, boundaries(k)));
}

void BinTree::Flip(IntegerVector pts){
  // Flips all point pt in pts from inactive to active or vice versa and adjusts loss.
  for(int pt:pts){
    IntegerVector pt_tree = orig_to_tree(_, pt);
    // for all i the following is the same (all active or all inactive). We choose i == 0
    bool type = tree(0, pt_tree(0));

    // Adjust loss part 1
    if(NodesOf(type) % 2 == 1)
      loss(type) -= sum(abs(FindMedian(type)-X(_, pt)));
    if(NodesOf(1-type) % 2 == 1)
      loss(1-type) += sum(abs(FindMedian(1-type)-X(_, pt)));

    // Remove pt from type and add it to 1-type (i.e. flip it)
    for(int l:Range(0,k))
      for(int i:Range(0,n-1))
        tree(i, boundaries(l) + floor(pt_tree(i)/pow(2, l))) += 1-2*type;


    // Adjust loss part 2 (note that odd before is even now and vice versa)
    if(NodesOf(type) % 2 == 1)
      loss(type) -= sum(abs(FindMedian(type)-X(_, pt)));
    if(NodesOf(1-type) % 2 == 1)
      loss[1-type] += sum(abs(FindMedian(1-type)-X(_, pt)));
  }
}

double BinTree::EvalStatL1(){
  // decrease of loss by splitting current setup of active and inactive nodes
  return(loss[2]-loss[1]-loss[0]);
}

NumericVector BinTree::FindMedian(bool type){
  // Returns median over all active (type=1) or inactive (type=0) points for each dimension
  int T_type = NodesOf(type);
  if(T_type==0){
    warning("no nodes of this type active");
    return(0);
  }
  NumericVector medians(n);
  if(T_type % 2 == 1){
    IntegerVector loc_mid = mthElement(type, (T_type-1)/2);
    for(int i:seqnn(0,n-1))
      medians(i) = (X(i, tree_to_orig(i, loc_mid(i))));
  }
  else{
    // Median is not unique
    // For the purpose of minimizing the loss it suffices to take any median
    IntegerVector loc_left = mthElement(type, T_type/2-1);
    // IntegerVector loc_right = mthElement(type, T_type/2);
    for(int i:seqnn(0,n-1))
      // medians(i) = (X(i, tree_to_orig(i, loc_left(i))) + X(i, tree_to_orig(i, loc_right(i))))/2;
      medians(i) = X(i, tree_to_orig(i, loc_left(i)));
  }
  return(medians);
}

IntegerVector BinTree::mthElement(bool type, int _m){
  // Uses tree structure to return m-th inactive (type=0) or active (type=1) element
  if(NodesOf(type)<=_m)
    stop("Not enough nodes of this type");
  IntegerVector positions(n);
  for(int i:seqnn(0,n-1)){
    int pos = 0;
    int m = _m;
    IntegerVector range = k-Range(1,k);
    for(int l:range){
      int left;
      if(type == 1)
        left = tree(i, boundaries[l]+2*pos);
      else
        left = pow(2,l) - tree(i, boundaries[l] + 2*pos);
      if(left > m)
        pos = 2*pos;
      else{
        pos = 2*pos+1;
        m -= left;
      }
    }
    positions(i) = pos;
  }

  return(positions);
}

List BinTree::OneStepSearchL1(int len){
  // Find maximal statistic of diagonal {(s,t): t-s=len}
  // returns s-1 and the maximal statistic
  // s-1 will be transformed to s later (see MaxStatsL1)
  // When calling [0:(len-1)] should be active (see BinTree::MaxStatsL1)
  int s = -1;
  double max = EvalStatL1();
  for(int i:seqnn(0,T-1-len)){
    Flip(IntegerVector::create(i, len+i));
    double stat = EvalStatL1();
    if(stat > max){
      max = stat;
      s = i;
    }
  }
  return(List::create(Named("shift") = s,
                      Named("stat") = max));
}

List BinTree::MaxStatsL1(String method="advanced", int min_seg=2){
  // method = "advanced" (optimistic search) or "full" (slow)
  // min_seg: minimal size of t-s
  // Returns the two best change point candidates ("shift", "ind")
  // and corresponding decrease in L1 loss "stat".
  // Note that to obtain the usual notion of "shift" and "ind",
  // one must add +1 (problem with C++ starting from 0). This is done in MaxStatsL1
  if(T<(2*min_seg))
    warning("ratio between maximal and minimal segment size to small");
  IntegerVector lengths;
  if(method=="advanced"){
    // lenghts are the diagonals that will be evaluated
    int lower = ceil(log2(min_seg));
    int upper = floor(log2(T-min_seg));
    lengths = vecpow(2, Range(lower,upper));
    if(lengths[0]!=min_seg)
      lengths.push_front(min_seg);
    if(*(lengths.end()-1) != T-min_seg)
      lengths.push_back(T-min_seg);
  }
  else if(method=="full"){
    // all diagonals will be evaluated
    lengths = seqnn(min_seg, T-min_seg);
  }else{warning("method must be either 'full' or 'advanced'");}
  int max_len = lengths[0];
  int max_shift = -1;
  double max_stat = 0;
  for(int len:lengths){
    ResetTree(true);
    Flip(Range(len_prev,len-1));
    len_prev = len;
    tree_prev = tree;
    loss_prev = loss;
    List cand = OneStepSearchL1(len);
    double stat = cand["stat"];
    if(stat>max_stat){
      max_stat = stat;
      max_len = len;
      max_shift = cand["shift"];
    }
  }
  int max_ind = max_shift + max_len;

// Linear search post-optimization for advanced search (new)
if(method=="advanced"){
  bool opt_shift = true;
  bool opt_ind = true;
  int iteration = 0;
  while((opt_shift || opt_ind) && iteration++<=3){
    if(opt_shift){
      // SHIFT OPTIMIZATION
      ResetTree(false);
      int shift_width = 1 + max_ind + 1 - min_seg;
      NumericVector shift_gains(shift_width);
      Flip(seqnn(max_ind-min_seg+1+1, max_ind));
      for(int i = max_ind-min_seg+1; i>= 0; i--){
        Flip(IntegerVector::create(i));
        shift_gains[i] = EvalStatL1();
      }
      int i = which_max(shift_gains);
      if(i != max_shift+1){
        max_shift = i-1;
        max_stat = shift_gains[i];
        opt_ind = true;
      }
      opt_shift = false;
    }

    if(opt_ind){
      // INDEX OPTIMIZATION
      ResetTree(false);
      int ind_width = T - (max_shift + min_seg);
      NumericVector ind_gains(ind_width);
      Flip(seqnn(max_shift+1, max_shift+min_seg-1));
      for(int i = max_shift+min_seg; i<T; i++){
        Flip(IntegerVector::create(i));
        ind_gains[i - (max_shift+min_seg)] = EvalStatL1();
      }
      int i = which_max(ind_gains);
      if(i + max_shift + min_seg != max_ind){
        max_ind = i + max_shift + min_seg;
        max_stat = ind_gains[i];
        opt_shift = true;
      }
      opt_ind = false;
    }
  }
}

  return(List::create(Named("shift") = max_shift,
                      Named("ind") = max_ind,
                      Named("stat") = max_stat));
}

List BinTree::BinarySegmentationL1(int min_seg = 2){
  // (non-circular) binary segmentation for L1 loss
  // Tree should be set to all nodes inactive i.e. BinTree::ResetTree(false);
  double max = 0;
  int ind=0;
  Flip(seqnn(0,min_seg-2));
  for(int i:Range(min_seg-1, T-1-min_seg)){
    Flip(IntegerVector::create(i));
    double stat = EvalStatL1();
    if(stat > max){
      max = stat;
      ind = i;
    }
  }
  return(List::create(Named("ind") = ind, Named("stat") = max));
}


List MaxStatsL1(NumericMatrix X, String method="advanced", int min_seg=2){
  // Make C++ code accessible from R. Adds +1 since integers in
  // C++ start from 0 but in R start from 1
  if(X.ncol()<min_seg)
    stop("not enough data");
  BinTree bt(X);
  List out = bt.MaxStatsL1(method, min_seg);
  int shift = out["shift"];
  int ind = out["ind"];
  return(List::create(Named("shift") = shift+1,
                      Named("ind") = ind+1,
                      Named("stat") = out["stat"]));
}

List BinarySegmentationL1(NumericMatrix X, int min_seg = 2){
  // Make C++ code accessible from R. Adds +1 since integers in
  // C++ start from 0 but in R start from 1
  BinTree bt(X);
  List out = bt.BinarySegmentationL1(min_seg);
  int ind = out["ind"];
  return(List::create(Named("shift") = 0, Named("ind") = ind+1, Named("stat") = out["stat"]));
}


double EvalStatL2(NumericMatrix S, int s, int t){
  // Returns gain from splitting the date at s and t.
  // "S" is cumulative sum matrix with S_i,0 = 0 for all dimensions i.
  const int T = S.ncol()-1;
  int n = S.nrow();
  if(s==t || (s==0 && t==T)){ // avoid div0
    warning("In EvalStatL2: s==t (div0)");
    return(0);
  }
  if(s>t){
    int temp = s;
    s = t;
    t = temp;
  }
  double res = 0;
  double inner;
  double outer;
  for(int i : Range(0,n-1)){
    inner = pow(S(i,t)-S(i,s),2)/(t-s);
    outer = pow(S(i,T)-S(i,t)+S(i,s),2)/(T-t+s);
    res += -pow(S(i,T),2)/T + inner  + outer;
  }
  if(res<0){
    // Rounding error; avoid crash from taking sqrt of negative number
    res=0;
  }
  return(sqrt(res));
}

// [[Rcpp::export]]
double EvalStatSlow(NumericMatrix X, int s, int t, String optimization="L2"){
  // Takes data matrix X as input. This function is slow and should be avoided
  // Value is the same as EvalStatL2 and BinTree::EvalStat
  // optimization = "L2":
  // Computes cumulative sums just for this evalutations
  // optimization = "L1":
  // Creates tree just for this evaluations
  int n = X.nrow();
  int T = X.ncol();
  if(optimization=="L2"){
    NumericMatrix S(n,T+1);
    for(int j:Range(0,n-1)){
      NumericVector Xj = X(j,_);
      NumericVector Sj = cumsum(Xj);
      Sj.push_front(0);
      S(j,_) = Sj;
    }
    return(EvalStatL2(S, s, t));
  }
  else if(optimization=="L1"){
    BinTree bt(X);
    bt.Flip(seqnn(s,t-1));
    return(bt.EvalStatL1());
  }
  stop("Optimization must be 'L1' or 'L2'");
  return 0;
}


List OneStepSearchL2(NumericMatrix S, int len){
  // Find maximal statistic of diagonal {(s,t): t-s=len}
  // returns s ("shift") and the maximal statistic ("stat")
  const int T = S.ncol()-1;
  double max_stat = EvalStatL2(S,0,len);
  int max_shift = 0;
  for(int i:seqnn(1, T-len)){
    double stat = EvalStatL2(S, i, i+len);
    if(stat > max_stat){
      max_stat = stat;
      max_shift = i;
    }
  }
  return(List::create(Named("shift") = max_shift, Named("stat") = max_stat));
}

List BinarySegmentationL2(NumericMatrix X, int min_seg = 2){
  // (non-circular) binary segmentation for L2 loss

  const int T = X.ncol();
  int n = X.nrow();

  NumericMatrix S(n, T+1);
  for(int j:Range(0,n-1)){
    NumericVector Xj = X(j,_);
    NumericVector Sj = cumsum(Xj);
    Sj.push_front(0);
    S(j,_) = Sj;
  }
  double max = 0;
  int ind=0;
  for(int i:seqnn(min_seg, T-min_seg)){
    double stat = EvalStatL2(S, 0, i);
    if(stat > max){
      max = stat;
      ind = i;
    }
  }
  return(List::create(Named("shift") = 0, Named("ind") = ind, Named("stat") = max));
}


List MaxStatsL2(NumericMatrix X, String method = "advanced",
                int min_seg = 2){
  // method = "advanced" (optimistic search) or "full" (slow)
  // min_seg: minimal size of t-s
  // Returns the two best change point candidates ("shift", "ind")
  // and corresponding decrease in L1 loss "stat".

  const int T = X.ncol();
  int n = X.nrow();

  NumericMatrix S(n, T+1);
  for(int j:Range(0,n-1)){
    NumericVector Xj = X(j,_);
    NumericVector Sj = cumsum(Xj);
    Sj.push_front(0);
    S(j,_) = Sj;
  }
  int max_seg = T-min_seg; // if len = T we evaluate the stat at (0,T) which is 0/0
  if(2*min_seg>T)
    warning("ratio between maximal and minimal segment size too small");
  IntegerVector lengths;
  if(method=="advanced"){
    int lower = ceil(log2(min_seg));
    int upper = floor(log2(max_seg));
    lengths = vecpow(2, Range(lower,upper));
    if(lengths(0)!=min_seg)
      lengths.push_front(min_seg);
    if(*(lengths.end()-1) != max_seg)
      lengths.push_back(max_seg);
  }
  else if(method=="full"){
    lengths = seqnn(min_seg, max_seg);
  }else{warning("method must be one of 'full' or 'advanced'");}
  int max_len = lengths(0);
  int max_shift = 0;
  double max_stat = 0;
  for(int len:lengths){
    List cand = OneStepSearchL2(S, len);
    double stat = cand["stat"];
    if(stat>max_stat){
      max_stat = stat;
      max_len = len;
      max_shift = cand["shift"];
    }
  }
  int max_ind = max_shift + max_len;

  // Linear search post-optimization
  if(method=="advanced"){
    bool shift_change = true;
    bool ind_change = true;
    int iteration = 0;
    while((shift_change || ind_change) && iteration++<=3 ){
      int search_range = (max_ind - max_shift)/2 +1;
      // SHIFT OPTIMIZATION
      // Boundaries: length of segment still >= 2, not under 1
      int shift_width_top = min(max_ind - max_shift - min_seg, search_range);
      int shift_width_bottom = min(min(max_shift, max_seg - (max_ind-max_shift)), search_range);
      int shift_width = shift_width_top + shift_width_bottom + 1;
      NumericVector shift_gains(shift_width);
      for(int i:seqnn(0, shift_width - 1))
        shift_gains(i) = EvalStatL2(S, max_shift - shift_width_bottom + i, max_ind);
      int i = which_max(shift_gains);
      if(i == shift_width_bottom)
        shift_change = false;
      else{
        max_shift += i - shift_width_bottom;
        ind_change = true;
      }

      search_range = (max_ind - max_shift)/2;
      // INDEX OPTIMIZATION
      // Boundaries: not over T, length of segment still >= 2
      int ind_width_top = min(min(T - max_ind, max_seg-(max_ind-max_shift)), search_range);
      int ind_width_bottom = min(max_ind - max_shift - min_seg, search_range);
      int ind_width = ind_width_top + ind_width_bottom + 1;
      NumericVector ind_gains(ind_width);

      for(int i:seqnn(0, ind_width - 1))
        ind_gains(i) = EvalStatL2(S, max_shift, max_ind - ind_width_bottom + i);
      i = which_max(ind_gains);
      if( i == ind_width_bottom)
        ind_change = false;
      else{
        max_ind += i - ind_width_bottom;

        shift_change = true;
      }
    }
    max_stat = EvalStatL2(S, max_shift, max_ind);
  }

  return(List::create(Named("shift") = max_shift,
                      Named("ind") = max_ind,
                      Named("stat") = max_stat));
}


//' Leading NA
//'
//' This function returns a logical vector identifying if
//' there are leading NA, marking the leadings NA as TRUE and
//' everything else as FALSE.
//'
//' @param x An integer vector
//' @export
// [[Rcpp::export]]
List MaxStats(NumericMatrix X, String optimization="L2", String method="advanced",
              bool circular=true, int min_seg=2){
  // Returns best change point candidates "shift","ind" and the corresponding
  // gain "stat".
  // For BS, the (interesting) change point is "ind" (and shift = 0)
  // optimization = "L1" or "L2", method = "advanced" or "full".
  // circular = 1 for CBS, circular = 0 for BS
  if(circular){
    if(optimization=="L1")
      return(MaxStatsL1(X, method, min_seg));
    else if(optimization=="L2")
      return(MaxStatsL2(X, method, min_seg));
  }
  else{
    if(optimization=="L1")
      return(BinarySegmentationL1(X, min_seg));
    else if(optimization=="L2")
      return(BinarySegmentationL2(X, min_seg));
  }
  warning("choose between L1 and L2 optimization");
  return(List::create());
}

List ChangePoints2(NumericMatrix X,
                   String optimization = "L2",
                   String method = "advanced"){
  // Returns exactly 2 change points which are always valid (no model selection).
  // Only used for presenting the results in the paper
  List cand = MaxStats(X, optimization, method, true, 2);
  IntegerVector cps = IntegerVector::create(cand["shift"], cand["ind"]);
  return(List::create(Named("cps") = cps,
                      Named("stat") = cand["stat"]));
}


double pPerm(NumericMatrix X, IntegerVector boundary, double cand_stat, String optimization, String method, double alpha,
              int nr_perms = 10000, bool circular = true, int min_seg=2){
  // Returns p value based on permutation tests. Uses early stopping.
  const int T = X.ncol();
  const int n = X.nrow();
  int r = ceil(nr_perms * alpha);
  double counter = 0;

  // Ignore early stopping by changing stopping boundary
  // boundary = rep(nr_perms,nr_perms*alpha)

  double p;
  for(int i:seqnn(1,nr_perms)){
    NumericMatrix X_sample(n, T);
    NumericVector Xj(T);
    for(int j:seqnn(0,n-1)){
      Xj = X(j,_);
      Xj = sample(Xj,T,false);
      X_sample(j,_) = Xj;
    }
    double perm_stat = MaxStats(X_sample, optimization, method, circular, min_seg)["stat"];

    if(cand_stat <= perm_stat){
      counter++;
      if(counter == r){
        // Rcout << i << " rejected with stat = " << cand_stat << endl;
        p = counter/i;
        return(p);
      }
    }
    if(boundary[counter]<i+1){
      // Rcout << i << " accepted with stat = " << cand_stat << endl;
      p = counter/i;
      return(p);
    }
    // if(i % 500 == 0)
    //   Rcout << i << " ";
  }
  // Rcout << " accepted with stat = " << cand_stat << endl;
  p = counter/nr_perms;
  return(p);
}


// [[Rcpp::export]]
bool PermTest(NumericMatrix X, IntegerVector boundary, double cand_stat, String optimization, String method, double alpha,
                int nr_perms = 10000, bool circular = true, int min_seg=2){
  // Does a permutation test using the p-value from pPerm function.
  double p = pPerm(X, boundary, cand_stat, optimization, method, alpha, nr_perms, circular, min_seg);
  if(p<alpha)
    return(true);
  else
    return(false);
}

//' Leading NA
//'
//' This function returns a logical vector identifying if
//' there are leading NA, marking the leadings NA as TRUE and asdfdsaf
//' everything else as FALSE. asdfdsaf
//'
//' @param x An integer vector
//' @keywords internal
//' @export
//'
// [[Rcpp::export]]
IntegerVector StoppingBoundary(int nr_perms = 10000, double alpha = 0.01, double eta_star = 0.0005) {
  // Calculates the early stopping boundary
  // eta_star is a control parameter for the probability of early stopping errors
  int r =  ceil(nr_perms*alpha);
  IntegerVector bdry(r);

  for(int i=0; i<r; i++){
    NumericVector y = dhyper(seqnn(0,i), bdry[i], nr_perms - bdry[i], r);
    while(sum(y) > eta_star){
      bdry[Range(i,r-1)] = bdry[Range(i,r-1)] + 1;
      y = dhyper(seq(0,i), bdry[i], nr_perms - bdry[i], r);
    }
  }
  return(bdry);
}




