#include <iostream>
#include <vector>
#include <unordered_map>

using namespace std;

//adjacency list of the section
vector<vector<int>> adj;
//section type vector T
vector<int> T;
//fault information vector I
vector<int> I;
//fault detection vector D
vector<int> D;
//fault candidate set X
vector<int> X;
//approximation gain ΔX
int DeltaX;

// Function to take user input
void takeUserInput() {
    int n;
    cout << "Enter the number of sections: ";
    cin >> n;

    cout << "Enter the adjacency list for the network topology:" << endl;
    adj.resize(n);
    for (int i = 0; i < n; i++) {
        cout << "Section " << i << ": ";
        int m;
        cin >> m;
        adj[i].resize(m);
        for (int j = 0; j < m; j++) {
            cin >> adj[i][j];
        }
    }

    cout << "Enter the section type vector T:" << endl;
    T.resize(n);
    for (int i = 0; i < n; i++) {
        cin >> T[i];
    }

    cout << "Enter the fault information vector I:" << endl;
    I.resize(n);
    for (int i = 0; i < n; i++) {
        cin >> I[i];
    }

    cout << "Enter the fault detection vector D:" << endl;
    D.resize(n);
    for (int i = 0; i < n; i++) {
        cin >> D[i];
    }
}

// Function to print the network topology as a connected graph
void printConnectedGraph() {
    cout << "Connected Network Topology:" << endl;
    for (int i = 0; i < adj.size(); i++) {
        cout << i;
        for (int j : adj[i]) {
            cout << " - " << j;
        }
        cout << endl;
    }
}

// Function to print the sections separately
void printSectionsSeparately() {
    cout << "Sections Separately:" << endl;
    for (int i = 0; i < adj.size(); i++) {
        cout << "Section " << i << ": ";
        for (int j = 0; j < adj[i].size(); j++) {
            cout << adj[i][j] << " ";
        }
        cout << endl;
    }
}

void DACFL(int i, int j) {
    if (i > j) {
        return;
    }

    int k = i;
    int Delta = I[k] - T[k] * T[k] * D[k];
    while (T[k] != 0 && Delta > 0) {
        k++;
        Delta = I[k] - T[k] * T[k] * D[k];
    }

    // Base case
    if (k == j) {
        X.clear();
        DeltaX = 0;
        for (int l = i; l <= j; l++) {
            if (T[l] == 1) {
                X.push_back(l);
                DeltaX += I[l] - T[l] * T[l] * D[l];
            }
        }
    } else {
        // Divide
        vector<int> u, v;
        for (int l = i; l <= k; l++) {
            if (T[l] == 1) {
                u.push_back(l);
            }
        }
        for (int l = k + 1; l <= j; l++) {
            if (T[l] == 1) {
                v.push_back(l);
            }
        }
        cout << "Dividing: Sections " << i << " to " << k << " and Sections " << k + 1 << " to " << j << endl;

        // Conquer
        int DeltaU = 0, DeltaV = 0;
        DACFL(i, k);
        DACFL(k + 1, j);

        // Combine
        X.clear();
        DeltaX = 0;
        vector<int> Xi_t;
        int DeltaXi_t = 0;
        for (int l = i; l <= j; l++) {
            if (T[l] == 1) {
                X.push_back(l);
                DeltaX += I[l] - T[l] * T[l] * D[l];
            } else {
                Xi_t.push_back(l);
                DeltaXi_t += I[l] - T[l] * T[l] * D[l];
            }
        }

        int maxDelta = DeltaXi_t;
        vector<int> maxXi_t = Xi_t;

        if (DeltaU > maxDelta) {
            maxDelta = DeltaU;
            maxXi_t = u;
        }
        if (DeltaV > maxDelta) {
            maxDelta = DeltaV;
            maxXi_t = v;
        }

        int DeltaU_V = 0;
        vector<int> X_u_v;
        for (int l = 0; l < u.size(); l++) {
            for (int m = 0; m < v.size(); m++) {
                int DeltaUV = I[u[l]] + I[v[m]] - T[u[l]] * T[u[l]] * D[u[l]] - T[v[m]] * T[v[m]] * D[v[m]];
                if (DeltaUV > DeltaU_V) {
                    X_u_v.clear();
                    X_u_v.push_back(u[l]);
                    X_u_v.push_back(v[m]);
                    DeltaU_V = DeltaUV;
                }
            }
        }

        if (DeltaU_V > maxDelta) {
            X = X_u_v;
            DeltaX = DeltaU_V;
        } else {
            X = maxXi_t;
            DeltaX = maxDelta;
        }
    }
}

int main() {
    takeUserInput();
    DACFL(0, T.size() - 1);
    printConnectedGraph();
    printSectionsSeparately();
    cout << "Fault candidate set X: ";
    for (int i = 0; i < X.size(); i++) {
        cout << X[i] << " ";
    }
    cout << endl;
    cout << "Approximation gain DeltaX: " << DeltaX << endl;
    return 0;
}


