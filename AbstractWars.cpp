#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <queue>

using namespace std;


class Base {
public:
    int m_id;
    int m_x;
    int m_y;
    int m_owner;
    int m_size;
    int m_gr; // growth rate
    bool m_under_attack;
};


ostream & operator<<(ostream & os, const Base &b){
    os << "Base - id: " << b.m_id 
       << " x: " << b.m_x 
       << " y: " << b.m_y 
       << " owner: " << b.m_owner
       << " size: " << b.m_size
       << " gr: " << b.m_gr
       << " ua: " << b.m_under_attack;

    return os;
}


class AttackParameters {
public:
    int boarder_base;
    int max_arival_time;

    AttackParameters() { boarder_base = -1; max_arival_time = -1; }
    AttackParameters(int bi, int mat): boarder_base(bi), max_arival_time(mat) {}
};


ostream & operator<<(ostream & os, const AttackParameters &ap){
    os << "AttackParameters - boarder_base: " << ap.boarder_base
       << " max_arival_time: " << ap.max_arival_time;

    return os;
}


double distance(Base &b1, Base &b2) {

    double x1 = b1.m_x;
    double y1 = b1.m_y;

    double x2 = b2.m_x;
    double y2 = b2.m_y;

    double xd = (x1 - x2) * (x1 - x2);
    double yd = (y1 - y2) * (y1 - y2);

    double d = sqrt( xd + yd );

    return d;
}


template<class T> void print_vector(vector<T> &v) {
    for (int i = 0; i < v.size(); ++i)
        cerr << v[i] << endl;
}


void print_matrix(vector<vector<double>> &dm) {

    fprintf(stderr, "Printing matrix...\n");
    for(int i = 0; i < dm.size(); i++) {
        for(int j = 0; j < dm[i].size(); j++)
            fprintf(stderr, "%4.2f ", dm[i][j]);
        fprintf(stderr, "\n");
    }

}


void update_bases(vector<int> &bases, vector<Base> &m_bases, int &m_step) {

    int N = bases.size()/2;

    // Calculate the growth rate on the first iteration.
    if (m_step == 1) {

        for(int i = 0; i < N; i++) {

            int gr = bases[2*i + 1] - m_bases[i].m_size;

            m_bases[i].m_owner = bases[2*i];
            m_bases[i].m_size = bases[2*i + 1];
            m_bases[i].m_gr = gr;
            m_bases[i].m_under_attack = false;
        }

    // On the other iterations updata only the owners and sizes.
    // Growth rate is constant.
    } else {

        for(int i = 0; i < N; i++) {
            m_bases[i].m_owner = bases[2*i];
            m_bases[i].m_size = bases[2*i + 1];
            m_bases[i].m_under_attack = false;
        }
    } // else end
}


int estimate_size_in_t_steps(Base &b, int &t_steps) {

    int current_size = b.m_size;
    int gr = b.m_gr;

    int future_size = current_size;

    for(int i = 0; i < t_steps; i++)
        future_size = future_size + (gr + future_size/100);

    if (future_size > 1000)
        return 1000;
    else
        return future_size;

}

int troop_arival_time(int &origin_base_id,
                      int &destination_base_id,
                      vector<vector<double>> m_destination_matrix,
                      int &speed) {

    double d = m_destination_matrix[origin_base_id][destination_base_id];
    return ceil(d/speed);
}



void estimate_attack_feasibility_for_base(AttackParameters &ap,
                                          vector<Base*> &base_neighbourhood,
                                          vector<vector<int>> &m_arival_time,
                                          int &max_number_of_attacking_bases) {

    int N = base_neighbourhood.size();
    int id = base_neighbourhood[0]->m_id;

    int number_of_checked_ally_bases = 0;
    int gathered_troops = 0;
    int boarder_base = -1;
    int max_arival_time = numeric_limits<int>::max();

    for(int i = 1; i < N; i++) {

        if (base_neighbourhood[i]->m_owner != 0)
            continue;

        if (number_of_checked_ally_bases > max_number_of_attacking_bases)
            break;

        // How much troops can we send from this base?
        gathered_troops = gathered_troops + (base_neighbourhood[i]->m_size)/2;

        // How long will it take the troops to arive at the enemies base.
        int attacking_base_id = base_neighbourhood[i]->m_id;
        int troop_arival_time = m_arival_time[id][attacking_base_id];      

        // What will be the size of the enemies base at troop arival time?
        int future_base_size = estimate_size_in_t_steps(*base_neighbourhood[0], troop_arival_time);

        cerr << "troop_arival_time: " << troop_arival_time << " gathered_troops: " << gathered_troops << " future_base_size: " << future_base_size << endl;

        number_of_checked_ally_bases++;
        // If there are enough troops attack.
        if (future_base_size < gathered_troops) {

            boarder_base = i;
            max_arival_time = troop_arival_time;
            break;

        } else {
            // If there are not enough troops continue gathering troops.
            continue;
        }
    }

    ap.boarder_base = boarder_base;
    ap.max_arival_time = max_arival_time;

}


void estimate_attack_feasibility(vector<AttackParameters> &attack_parameters_per_base,
                                 vector<vector<Base*>> &base_neighbourhood,
                                 vector<vector<int>> &m_arival_time,
                                 int &max_number_of_attacking_bases) {

    int N = base_neighbourhood.size();

    for(int i = 0; i < N; i++) {
        cerr << "Processing base_neighbourhood i:" << i << endl;
        // When attacking consider only enemy bases.
        if (base_neighbourhood[i][0]->m_owner != 0)
            continue;

        estimate_attack_feasibility_for_base(attack_parameters_per_base[i],
                                             base_neighbourhood[i],
                                             m_arival_time,
                                             max_number_of_attacking_bases);

    }


}


class AbstractWars {
public:
    static int m_step;
    int m_speed;
    vector<Base> m_bases;
    vector<vector<double>> m_distance_matrix;
    vector<vector<int>> m_arival_time;
    vector<vector<Base*>> m_base_neighbourhood;
    vector<AttackParameters> m_attack_parameters_per_base;

    AbstractWars() {};
    int init(vector <int> base_locations, int speed);
    vector <int> sendTroops(vector <int> bases, vector <int> troops);
};


int AbstractWars::init(vector <int> base_locations, int speed) {

    int N = base_locations.size()/2;

    cerr << "Number of base locations: " << N << endl;

    m_bases.resize(N);
    m_distance_matrix.resize(N);
    m_arival_time.resize(N);
    m_base_neighbourhood.resize(N);
    m_attack_parameters_per_base.resize(N);
    m_speed = speed;

    for(int i = 0; i < N; i++) {
        m_bases[i].m_id = i;
        m_bases[i].m_x = base_locations[2*i];
        m_bases[i].m_y = base_locations[2*i + 1];
        m_bases[i].m_under_attack = false;
        m_distance_matrix[i].resize(N);
        m_arival_time[i].resize(N);
        m_base_neighbourhood[i].resize(N);
    }


    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {

            m_base_neighbourhood[i][j] = &m_bases[j];

            if (i <= j)
                continue;

            double d = distance(m_bases[i], m_bases[j]);
            m_distance_matrix[i][j] = d;
            m_distance_matrix[j][i] = d;

            int t = ceil(d/speed);
            m_arival_time[i][j] = t;
            m_arival_time[j][i] = t;

        }
    }



    for(int i = 0; i < m_bases.size(); i++) {

        auto cf = [this, &i](const Base *b1, const Base *b2){ 

                int id1 = b1->m_id;
                int id2 = b2->m_id;

                return m_distance_matrix[i][id1] < m_distance_matrix[i][id2]; 

            };

        sort(m_base_neighbourhood[i].begin(), m_base_neighbourhood[i].end(), cf);

    }

    cerr << "Print bases: " << endl;
    print_vector(m_bases);

    cerr << "Print distance matrix:" << endl;
    print_matrix(m_distance_matrix);

    return 0;
}


vector<int> AbstractWars::sendTroops(vector <int> bases, vector <int> troops) {

    cerr << "Step: " << m_step << endl;
    update_bases(bases, m_bases, m_step);
    print_vector(m_bases);


    int t_steps = 10;
    int index = 3;
    int fs = estimate_size_in_t_steps(m_bases[index], t_steps);

    cerr << "current: " << m_bases[index].m_size << " fs: " << fs << endl;

    int orig = 0;
    int dest = 10;

    cerr << "troop arival time: " << troop_arival_time(orig, dest, m_distance_matrix, m_speed) << endl;
    cerr << m_arival_time[orig][dest] << endl;
        
    index = 4;
    for(int j = 0; j < m_bases.size(); j++) {

        int id1 = m_base_neighbourhood[index][j]->m_id;
        cerr << *m_base_neighbourhood[index][j] << " d: " << m_distance_matrix[index][id1] << endl;
    
    }

    int max_number_of_attacking_bases = 100;

    estimate_attack_feasibility(m_attack_parameters_per_base,
                                m_base_neighbourhood,
                                m_arival_time,
                                max_number_of_attacking_bases);

    print_vector(m_attack_parameters_per_base);

    vector<int> res;



    m_step++;
    return res;
}

int AbstractWars::m_step = 0;

// -------8<------- end of solution submitted to the website -------8<-------

template<class T> void getVector(vector<T>& v) {
    for (int i = 0; i < v.size(); ++i)
        cin >> v[i];
}

int main() {
    AbstractWars aw;
    int N;
    cin >> N;
    vector<int> baseLoc(N);
    getVector(baseLoc);
    int S;
    cin >> S;
    
    int retInit = aw.init(baseLoc, S);
    cout << retInit << endl;
    cout.flush();
    
    for (int st = 0; st < 2000; ++st) {
        int B;
        cin >> B;
        vector<int> bases(B);
        getVector(bases);
        int T;
        cin >> T;
        vector<int> troops(T);
        getVector(troops);
        
        vector<int> ret = aw.sendTroops(bases, troops);
        cout << ret.size() << endl;
        for (int i = 0; i < (int)ret.size(); ++i) {
            cout << ret[i] << endl;
        }
        cout.flush();
    }
}
