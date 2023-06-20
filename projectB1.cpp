#include <bits/stdc++.h>
using namespace std;

class Graph
{
private:
    vector<int> affected; // list of affected nodes
    vector<set<int>> adjList;
    vector<int> oldSize;    // r(y)
    vector<int> newSize;    // r'(y)
    vector<int> compNumber; // to calculate if u and v are in same component
    vector<float> closeness;
    vector<float> estimated;
    vector<int> cutoffDist;
    vector<bool> isExact;
    int numVertices;
    vector<bool> changedReach;
    set<pair<float, int>> TopK;
    set<pair<float, int>> deleteQueue;
    int numEdges;
    int diameter;
    int k;

public:
    Graph()
    {
        numVertices = 1000;
        numEdges = 0;
        diameter = 0;
        k = 20;
        adjList.resize(1000);
        closeness.resize(1000, 0);
        estimated.resize(1000, 0);
        cutoffDist.resize(1000, 0);
        isExact.resize(1000, false);
        changedReach.resize(1000, false);
        oldSize.resize(1000, 0);
        newSize.resize(1000, 0);
        compNumber.resize(1000, 0);
        for (int i = 0; i < numVertices; i++)
        {
            oldSize[i] =1; newSize[i] = 1;
            compNumber[i] = i;
        }
        
    }
    //copy constructor
    Graph(const Graph &G)
    {
        numVertices = G.numVertices;
        numEdges = G.numEdges;
        diameter = G.diameter;
        k = G.k;
        adjList = G.adjList;
        oldSize = G.oldSize;
        newSize = G.newSize;
        compNumber = G.compNumber;
        closeness = G.closeness;
        estimated = G.estimated;
        cutoffDist = G.cutoffDist;
        isExact = G.isExact;
        changedReach = G.changedReach;
        TopK = G.TopK;
        deleteQueue = G.deleteQueue;
    }

    vector<int> BFS(int root) // BFS returning the PM
    {
        // normal BFS to calculate the distance
        vector<bool> visited(numVertices, false);
        vector<int> distance(numVertices, 1e8);
        queue<int> Q;
        distance[root] = 0;
        Q.push(root);
        while (!Q.empty())
        {
            int vertex = Q.front();
            Q.pop();
            if (visited[vertex])
                continue;
            visited[vertex] = true;

            for (auto child : adjList[vertex])
            {
                if (visited[child])
                    continue;
                Q.push(child);
                if (distance[child] > distance[vertex] + 1)
                    distance[child] = distance[vertex] + 1;
            }
        }

        return distance;
    }

    vector<float> BFSCut(int node, float xk) // To be implemented
    {
        float ccUpperbound = adjList[node].size() + (newSize[node] - adjList[node].size()) / 2.0;
        vector<float> ans(3);
        if (ccUpperbound < xk)
        {
            ans[0] = ccUpperbound;
            ans[1] = 0;
            ans[2] = 0;
            return ans;
        }
        vector<bool> visited(numVertices, false);
        vector<int> distance(numVertices, 1e8);
        queue<int> Q;
        distance[node] = 0;
        int nodesvisited = 0;
        Q.push(node);

        int currentdepth = 0;
        float harmonicSum = 0;
        float extraHarmonicSum = 0;
        int estimatedNodesatDepthPlus1 = 0;
        while (!Q.empty())
        {
            int vertex = Q.front();
            Q.pop();
            if (visited[vertex])
                continue;
            visited[vertex] = true;
            nodesvisited++;

            for (auto child : adjList[vertex])
            {
                if (visited[child])
                    continue;
                Q.push(child);
                if (distance[child] > distance[vertex] + 1)
                    distance[child] = distance[vertex] + 1;
                if (distance[child] == currentdepth)
                {
                    extraHarmonicSum += 1.0 / distance[child];
                    estimatedNodesatDepthPlus1 += adjList[child].size();
                    nodesvisited++;
                }
                else if (distance[child] == currentdepth + 1)
                {
                    float newcc = (harmonicSum + extraHarmonicSum) + (estimatedNodesatDepthPlus1) / (currentdepth + 1.0) + (newSize[node] - estimatedNodesatDepthPlus1 - nodesvisited) / (currentdepth + 2.0);
                    if (newcc < xk)
                    {
                        ans[0] = newcc;
                        ans[1] = 0;
                        ans[2] = currentdepth;
                        return ans;
                    }
                    else
                    {
                        currentdepth++;
                        harmonicSum += extraHarmonicSum;
                        extraHarmonicSum = 1.0 / distance[child];
                        estimatedNodesatDepthPlus1 = adjList[child].size();
                        nodesvisited++;
                    }
                }
            }
        }
        float newcc = (harmonicSum + extraHarmonicSum) + (estimatedNodesatDepthPlus1) / (currentdepth + 1.0) + (newSize[node] - estimatedNodesatDepthPlus1 - nodesvisited) / (currentdepth + 2.0);
        ans[0] = newcc;
        ans[1] = 1;
        ans[2] = currentdepth;
        return ans;
    }

    pair<vector<int>, vector<int>> findAffected(int v1, int v2, string instruction)
    {
        // clearing the orignal affected list
        affected.clear();

        // calculating the distances before the addition
        vector<int> oldDistv1 = BFS(v1);
        vector<int> oldDistv2 = BFS(v2);

        // adding the edge
        if (instruction == "add")
        {
            adjList[v1].insert(v2);
            adjList[v2].insert(v1);
        }
        if (instruction == "delete")
        {
            adjList[v1].erase(v2);
            adjList[v2].erase(v1);
        }

        // calculating the distance after the addition (reach is also updated here)
        vector<int> newDistv1 = BFS(v1);
        vector<int> newDistv2 = BFS(v2);

        // if distances are changed then they are affected
        for (int i = 0; i < numVertices; i++)
        {
            if (oldDistv1[i] != newDistv1[i] || oldDistv2[i] != newDistv2[i])
            {
                affected.push_back(i);
            }
        }

        return {newDistv1, newDistv2};
    }

    void BFSMark(int root, int mark) // mark the nodes into a component with compNumber = mark
    {
        vector<bool> visited(numVertices, false);
        queue<int> Q;
        Q.push(root);
        while (!Q.empty())
        {
            int vertex = Q.front();
            Q.pop();
            if (visited[vertex])
                continue;
            visited[vertex] = true;
            compNumber[vertex] = mark;
            for (auto child : adjList[vertex])
            {
                if (visited[child])
                    continue;
                Q.push(child);
            }
        }
    }

    void updateReach(int v1, int v2, string instruction) // updating the Reach
    {
        if (instruction == "add")
        {
            // if they are in different components
            if (compNumber[v1] != compNumber[v2])
            {
                // shifting component 2 to component 1

                // updating sizes
                oldSize[v1] = newSize[v1];
                oldSize[v2] = newSize[v2];
                newSize[v1] = oldSize[v1] + oldSize[v2];
                newSize[v2] = 0;
                int changeFrom = compNumber[v2];
                int changeTo = compNumber[v1];

                // shifting the vertices
                for (int i = 0; i < numVertices; i++)
                {
                    if (compNumber[i] == changeFrom)
                    {
                        compNumber[i] = changeTo;
                    }
                }
            }
        }
        if (instruction == "delete")
        {
            // assuming the edge is deleted
            BFSMark(v1, v1); // marking the nodes in the reach of v1 as v1
            BFSMark(v2, v2); // marking the nodes in the reach of v2 as v2
        }
    }

    void processAffected(int node, float &kSmallest, int u, vector<int> &distU, vector<int> &oldDepthSize, vector<int> &newDepthSize)
    {
        // vertex node is affected and we have the smallest value in Top K list

        // condition 1
        bool isFar = (cutoffDist[node] < distU[node]);
        // condition 2
        bool onBorder = (cutoffDist[node] == distU[node]);

        if (isFar && !isExact[node]) // condition 1
        {
            estimated[node] += (1.0 * (newSize[node] - oldSize[node])) / (1.0 * (2 + cutoffDist[node]));
        }
        else if (onBorder && !isExact[node])
        {
            estimated[node] += (1.0 * (cutoffDist[node] + 1) + (1.0 * (newSize[node] - oldSize[node] - 1)) / (1.0 * (2 + cutoffDist[node])));
        }
        else
        {
            // updateClosenessTotal(node);
            float additionalValue = 0;
            for (int i = 1; i <= newDepthSize.size(); i++)
            {
                additionalValue += ((1 / (i + distU[node])) * newDepthSize[i]);
            }
            for (int i = 1; i <= oldDepthSize.size(); i++)
            {
                additionalValue -= ((1 / (i + distU[node])) * oldDepthSize[i]);
            }
            cutoffDist[node] = newDepthSize.size();
        }

        // checking if we need to process further
        if (estimated[node] >= kSmallest)
        {
            // Using the BFS Cut for the node
            vector<float> retBuffer = BFSCut(node, kSmallest);
            estimated[node] = retBuffer[0];
            isExact[node] = retBuffer[1];
            cutoffDist[node] = retBuffer[2];

            // checking on the updated parameters if we can add it to the TopK list
            if (isExact[node] && estimated[node] > kSmallest)
            {
                // adding to the list
                TopK.insert({estimated[node], node});

                // if size exceeds then remove the least one
                if (TopK.size() > k)
                {
                    TopK.erase(TopK.begin());
                }

                // if size is equal to k then we update the smallest value here
                if (TopK.size() == k)
                {
                    kSmallest = (*TopK.begin()).first;
                }
            }
        }
    }

    void updateClosenessTotal(int node) // updating the closeness of a given node (Exact value)
    {
        vector<int> distance = BFS(node);
        double value = 0; // the closeness
        for (int j = 0; j < numVertices; j++)
        {
            if (j == node || distance[j] == 1e8)
                continue;

            value += (1.0 / (1.0 * distance[j]));
        }
        closeness[node] = value;
        isExact[node] = true;
        cutoffDist[node] = diameter;
    }

    void cleanTopKList()
    {
        // if a affected node is present in the TopK list then remove it
        for (auto node : affected)
        {
            if (TopK.find({estimated[node], node}) == TopK.end())
                continue;
            TopK.erase({estimated[node], node});
        }
    }

    vector<int> findNumDepth(int node) // returns a vector which stores the number of nodes at distance i
    {
        vector<int> distance = BFS(node);
        int maxDepth = 0;
        for (int i = 0; i < numVertices; i++)
        {
            if (distance[i] == 1e8)
                continue;
            maxDepth = max(maxDepth, distance[i]);
        }

        vector<int> numDepths(maxDepth + 1, 0);
        for (int i = 0; i < numVertices; i++)
        {
            if (distance[i] == 1e8)
                continue;
            numDepths[distance[i]]++;
        }

        return numDepths;
    }

    void fillDeleteQueue() // the Queue used in the delete edge operation is filled
    {
        for (int i = 0; i < numVertices; i++)
        {
            deleteQueue.insert({estimated[i], i});
        }
    }

    float getTopK(int index) // returns the index th index in the set
    {
        auto it = TopK.begin();
        for (int i = 0; i < index; i++)
        {
            it++;
        }

        return (*it).first;
    }

    void processQueue(float &kSmallest) // processes the queue used for deletion
    {
        while (!deleteQueue.empty())
        {
            pair<float, int> information = *deleteQueue.rbegin(); // the max element
            deleteQueue.erase(information);                       // popping the max element
            int node = information.second;
            int value = information.first;

            // when size in more
            if (TopK.size() >= k)
            {
                float smallK = getTopK(TopK.size() - k); // TopK[K] -> kth smallest
                if (value > smallK)
                {
                    return;
                }
            }

            // if its not the exact node
            if (!isExact[node])
            {
                // Using the BFS Cut for the node
                vector<float> retBuffer = BFSCut(node, kSmallest);
                estimated[node] = retBuffer[0];
                isExact[node] = retBuffer[1];
                cutoffDist[node] = retBuffer[2];
            }

            // if its the exact node and the value is more than the k smallest
            if (isExact[node] && value > kSmallest)
            {
                TopK.insert(information);
                if (TopK.size() > k)
                {
                    TopK.erase(TopK.begin());
                }
                if (TopK.size() == k)
                {
                    kSmallest = (*TopK.begin()).first;
                }
            }
        }
    }

    void addEdge(int v1, int v2) // adding an undirected edge between v1 and v2
    {
        // finding the old sizes of the depth list
        vector<int> numOldDepthv1 = findNumDepth(v1);
        vector<int> numOldDepthv2 = findNumDepth(v2);

        // finding the affected list and getting the new Distances
        pair<vector<int>, vector<int>> newDistances;
        newDistances = findAffected(v1, v2, "add");
        vector<int> newDistv1 = newDistances.first;
        vector<int> newDistv2 = newDistances.second;

        // updating the reach
        updateReach(v1, v2, "add");

        // removing affected nodes from TopK list
        cleanTopKList();

        // get the least values node in the list
        float xk = (*TopK.begin()).first;

        // finding the new sizes of the depth list
        vector<int> numNewDepthv1 = findNumDepth(v1);
        vector<int> numNewDepthv2 = findNumDepth(v2);

        // process the affected nodes
        for (auto node : affected)
        {
            // Need to update reach, estimation, and value of ni
            processAffected(node, xk, v1, newDistv1, numOldDepthv1, numNewDepthv1);
            processAffected(node, xk, v2, newDistv2, numOldDepthv2, numNewDepthv2);
        }
    }

    void deleteEdge(int v1, int v2)
    {
        // finding the affected list and getting the new Distances
        pair<vector<int>, vector<int>> newDistances;
        newDistances = findAffected(v1, v2, "delete");
        vector<int> newDistv1 = newDistances.first;
        vector<int> newDistv2 = newDistances.second;

        // updating the reach
        updateReach(v1, v2, "delete");

        // removing affected nodes from TopK list
        cleanTopKList();

        // initializing the kth smallest number
        float kSmallest = 0;

        // Filling the priority Queue (set here) for deletion algo
        fillDeleteQueue();

        processQueue(kSmallest);
    }

    void printGraph()
    {
        for (int i = 0; i < numVertices; i++)
        {
            cout << i << " : ";
            for (auto child : adjList[i])
            {
                cout << child << " ";
            }
            cout << endl;
        }
    }
};

int main()
{
    Graph G;
    G.addEdge(0, 1);
    G.addEdge(1, 2);
    G.printGraph();
}