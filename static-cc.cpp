#include <bits/stdc++.h>

using namespace std;

class Graph{
    private:
        int numVertices;
        vector<vector<int>> E;
        int k = 5;
    public:
        Graph(int n){
            numVertices=n;
            E.resize(n);
        }
        void AddEdge(int u,int v){
            E[u].push_back(v);
            E[v].push_back(u);
        }
        vector<float> CC(){
            vector<float> cc;
            for(int i=0;i<numVertices;i++){
                map<int,bool> mark;
                int dist[numVertices];
                dist[i]=0;
                int sum=0;
                queue<int> q;
                q.push(i);
                while(!q.empty()){
                    int curr=q.front();
                    q.pop();
                    if(!mark[curr]){
                        sum+=dist[curr];
                        mark[curr]=true;
                        for(int i=0;i<E[curr].size();i++){
                            q.push(E[curr][i]);
                            dist[E[curr][i]]=dist[curr]+1;
                        }
                    }
                }
                if(sum==0) cc.push_back(0); //in case of isolated vertices
                else cc.push_back(1/(float)sum);
            }
            return cc;
        }

        vector<pair<int,float>> find_top_k() {

          vector<float> cc = CC();
          vector<pair<float,int>> sorted_cc;
          vector<pair<int,float>> top_k;

          for ( int i = 0 ; i < cc.size() ; i++ ) {

            sorted_cc.push_back({cc[i],i});

          }

          sort(sorted_cc.begin(),sorted_cc.end(),greater<>());

          for(int i = 0 ; i < k ; i++ ) {

            top_k.push_back({sorted_cc[i].second,sorted_cc[i].first});

          }

          return top_k;

        }

};

int main() {
    int n,m;
    cin>>n>>m;
    Graph G=Graph(n);
    for(int i=0;i<m;i++){
        int u,v;
        cin>>u>>v;
        G.AddEdge(u,v);
    }
    vector<float> v=G.CC();
    for(int i=0;i<v.size();i++){
        cout<<v[i]<<" ";
    }
}
