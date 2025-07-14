#include<vector>
#include<random>
#include<iostream>
#include<algorithm>
#include<cmath>
using namespace std;
typedef long long int ll;

int randomInt(int n) {      // Returns a integer in {0,...,n} generated uniformly at random
    static random_device rd; static mt19937 gen(rd());
    uniform_int_distribution<> dist(0, n);
    return dist(gen);
}

class Tree {
    protected:
        int size = 0;       // Number of vertices in the tree
        vector<int> parent = {};     // i'th slot of this vector contains the index of the parent of the i'th added vertex (root has no parent so we put -1)
        vector<vector<int>> children = {};     // i'th slot of this vector contains a list of the children of the i'th added vertex
        bool phiToDate = true;      // indicates whether the phi vector is up to date
        vector<double> logPhi = {};      // (if phiToDate = true) i'th slot of this vector contains log(phi(i))
    public:
        Tree() = default;
        Tree(vector<int> parents) {       // You can use this constructor if you want your tree to have an exact structure. logPhi will not be automatically updated!
            size = parents.size();
            parent = parents;
            children = vector<vector<int>>(size);
            for (int i = size-1; i > 0; i--) children[parent[i]].push_back(i);
            phiToDate = false;
            logPhi = vector<double>(size);
        }
        void add(int par) {      // Adds a vertex to the tree that will have par as its parent
            parent.push_back(par);
            children.push_back({});
            children[par].push_back(size++);
            phiToDate = false;
        }
        void updatePhi() {      // O(n) time algorithm for computing all values of Phi, only uses one vector of length n of buffer memory
            if (size == 0) return;
            phiToDate = true;
            vector<double> &sizeSubtree = logPhi;        // i'th slot of this vector contains the number of descendants of vertex i (including i) w.r.t. the root
            for (int i = 0; i < size; i++) sizeSubtree[i] = 1;
            for (int i = size-1; i > 0; i--) {      // sizeSubtree[i] = 1 + \sum_{j is a child of i} sizeSubtree[j]
                sizeSubtree[parent[i]] += sizeSubtree[i];
            }

            vector<double> logProdSubtree;        // i'th slot of this vector contains the product of the sizes of all subtrees of descendants of vertex i (including i) w.r.t. the root
            for (int k : sizeSubtree) logProdSubtree.push_back(log(k));
            for (int i = size-1; i > 0; i--) {       // logProdSubtree[i] = log(sizeSubtree[i]) + \sum_{j is a child of i} logProdSubtree[j]
                logProdSubtree[parent[i]] += logProdSubtree[i];
            }

            vector<double> &logProdSuptree = sizeSubtree;       // i'th slot of this vector contains the product of the sizes of all subtrees needed to calculate phi(i) that aren't rooted at a descendant of i w.r.t. to root
            logProdSuptree[0] = 0.0;
            for (int i = 1; i < size; i++) {        // logProdSuptree[i] = log(size - sizeSubtree[i]) + logProdSuptree[parent[i]] - logProdSubtree[i] + \sum_{j is a child of parent[i]} logProdSubtree[j]
                logProdSuptree[i] =  log(size - sizeSubtree[i]) + logProdSuptree[parent[i]] - logProdSubtree[i];
                for (int j : children[parent[i]]) logProdSuptree[i] += logProdSubtree[j];
            }

            for (int i = 0; i < size; i++) for (int j : children[i]) logPhi[i] += logProdSubtree[j];       // logPhi[i] = logProdSuptree[i] + \sum_{j is a child of i} logProdSubtree[j]
        }
        vector<double> getLogPhi() {
            return logPhi;
        }
        void printLogPhi() {
            cout << "{";
            if (!(logPhi.empty())) {
                cout << logPhi[0];
                for (int i = 1; i < size; i++) cout << ", " << logPhi[i];
            }
            cout << "}" << endl;
        }
        bool rootInTop(int k) {       // Returns true if and only if there are at most k-1 vertices i other than 0 such that phi(i) <= phi(0)
            if (k >= size) return true;
            vector<double> tmp = logPhi;
            sort(tmp.begin(),tmp.end());
            if (logPhi[0] < tmp[k]) return true;
            return false;
        }
};

class UATree : public Tree {
    public:
        UATree(int n) {     // Builds a tree of size n by uniform attachment. logPhi will not be automatically updated!
            for (int i = 0; i < n; i++) add();
        }
        void add() {       // Adds a vertex by uniform attachment
            if (size == 0) parent.push_back(-1);
            else {
                int par = randomInt(size-1);
                parent.push_back(par);
                children[par].push_back(size);
            }
            size++;
            children.push_back({});
            phiToDate = false;
            logPhi.push_back(0.0);
        }
};

ll K(double c, double eps) {       // This is the function K from the paper
    return ceil(exp(c + c*sqrt(log(1.0/eps))));
}

double prob(ll K, int n, int N) {       // Outputs an approximation of the probability that the root is among the top K lowest phi value vertices of a UA(n) distributed tree. The bigger N is, the better the approximation, but the more time it takes
    int cnt = 0;
    for (int i = 0; i < N; i++) {
        UATree tree = UATree(n);
        tree.updatePhi();
        if (tree.rootInTop(K)) cnt++;
    }
    return ((double) cnt) / ((double) N);
}

bool solve(double c, int n) {       // Checks (approximately) whether the theorem of our paper holds with c is the constant
    double eps = 1.0;
    while (eps > 2e-4) {
        int Ka = K(c,eps);
        if (Ka < 0) Ka = 2147483647;
        if (prob(Ka, n, 25/eps) < 1 - eps) return false;
        eps /= 4;
    }
    return true;
}

int main() {
    for (int i = 0; i < 20; i++) {      // Doing this 20 times
        double lo = 0.0;
        double hi = 3.0;
        while (!(solve(hi,1000))) {hi *= 2; lo = hi/2;}     // Finding a valid value for hi

        double mid = (lo + hi)/2;
        while (hi - lo > 1e-3) {       // Binary searching for c
            mid = (lo + hi)/2;
            cout << "mid = " << mid << endl;
            if (solve(mid, 1000)) hi = mid;
            else lo = mid;
        }
        cout << "c = " << (lo + hi)/2 << endl;
    }
}