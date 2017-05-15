#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <queue>



using namespace std;

const int LOCAL_MAX_INT = numeric_limits<int>::max();


class Base {
public:
    int m_id;
    int m_x;
    int m_y;
    int m_owner;
    int m_size;
    int m_gr; // growth rate
    bool m_under_attack;
    bool m_attacking;
    int m_attack_time;
};


ostream & operator<<(ostream & os, const Base &b){
    os << "Base - id: " << b.m_id 
       << " x: " << b.m_x 
       << " y: " << b.m_y 
       << " owner: " << b.m_owner
       << " size: " << b.m_size
       << " gr: " << b.m_gr
       << " ua: " << b.m_under_attack
       << " at: " << b.m_attack_time;

    return os;
}


class AttackParameters {
public:
    int m_destination_base_id;
    int m_destination_base_gr;
    int m_boarder_base;
    int m_max_arival_time;

    AttackParameters() { m_destination_base_id = -1; m_destination_base_gr = -1; m_boarder_base = -1; m_max_arival_time = -1; }

    AttackParameters(int dest, int gr, int bi, int mat): m_destination_base_id(dest), 
                                                         m_destination_base_gr(gr),
                                                         m_boarder_base(bi), 
                                                         m_max_arival_time(mat) {}

    AttackParameters(const AttackParameters &ap): m_destination_base_id(ap.m_destination_base_id),
                                                  m_destination_base_gr(ap.m_destination_base_gr), 
                                                  m_boarder_base(ap.m_boarder_base), 
                                                  m_max_arival_time(ap.m_max_arival_time) { 
                                                    //cerr << "Calling AttackParameters copy constructor." << endl;
                                                  }
};


ostream & operator<<(ostream & os, const AttackParameters &ap){
    os << "AttackParameters - m_destination_base_id: " << ap.m_destination_base_id
       << " m_destination_base_gr: " << ap.m_destination_base_gr
       << " boarder_base: " << ap.m_boarder_base
       << " max_arival_time: " << ap.m_max_arival_time;

    return os;
}


class PriorityQueueAttackParametersCompare {
	public:
		inline bool operator ()(const AttackParameters &ap1, const AttackParameters &ap2) const {
			return ap1.m_destination_base_gr < ap2.m_destination_base_gr;
		}
};


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
    for (int i = 0; i < v.size(); ++i) {
        cerr << v[i] << endl;
    }
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
            m_bases[i].m_attacking = false;
            m_bases[i].m_attack_time = -1;
        }

    // On the other iterations updata only the owners and sizes.
    // Growth rate is constant.
    } else {

        for(int i = 0; i < N; i++) {

            if ( (m_bases[i].m_under_attack == true) && (m_bases[i].m_attack_time == m_step) ) {
                m_bases[i].m_under_attack = false;
                m_bases[i].m_attack_time = -1;
            }

            m_bases[i].m_owner = bases[2*i];
            m_bases[i].m_size = bases[2*i + 1];
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

int troop_arrival_time(int &origin_base_id,
                      int &destination_base_id,
                      vector<vector<double>> m_destination_matrix,
                      int &speed) {

    double d = m_destination_matrix[origin_base_id][destination_base_id];
    return ceil(d/speed);
}



AttackParameters estimate_attack_feasibility_for_base(vector<Base*> &base_neighbourhood,
                                                      vector<vector<int>> &m_arival_time,
                                                      int &max_number_of_attacking_bases) {

    int N = base_neighbourhood.size();
    int id = base_neighbourhood[0]->m_id;
    int gr = base_neighbourhood[0]->m_gr;

    int number_of_checked_ally_bases = 0;
    int gathered_troops = 0;
    int boarder_base = -1;
    int max_arival_time = LOCAL_MAX_INT;

    for(int i = 1; i < N; i++) {

        if (base_neighbourhood[i]->m_owner != 0)
            continue;

        if (number_of_checked_ally_bases > max_number_of_attacking_bases)
            break;

        // How much troops can we send from this base?
        gathered_troops = gathered_troops + (base_neighbourhood[i]->m_size)/2;

        // How long will it take the troops to arive at the enemies base.
        int attacking_base_id = base_neighbourhood[i]->m_id;
        int troop_arrival_time = m_arival_time[id][attacking_base_id];      

        // What will be the size of the enemies base at troop arival time?
        int future_base_size = estimate_size_in_t_steps(*base_neighbourhood[0], troop_arrival_time);

        //cerr << "troop_arrival_time: " << troop_arrival_time << " gathered_troops: " << gathered_troops << " future_base_size: " << future_base_size << " growth rate: " << base_neighbourhood[0]->m_gr << endl;

        number_of_checked_ally_bases++;
        // If there are enough troops attack.
        if (1.2*future_base_size < gathered_troops) {

            boarder_base = i;
            max_arival_time = troop_arrival_time;
            break;

        } else {
            // If there are not enough troops continue gathering troops.
            continue;
        }
    }

    return AttackParameters(id, gr, boarder_base, max_arival_time);
}


AttackParameters estimate_attack_feasibility(vector<vector<Base*>> &base_neighbourhood,
                                             vector<vector<int>> &m_arival_time,
                                             int &max_number_of_attacking_bases,
                                             vector<Base> &bases) {

    int N = base_neighbourhood.size();

    priority_queue<AttackParameters, vector<AttackParameters>, PriorityQueueAttackParametersCompare> pq;

    for(int i = 0; i < N; i++) {

        //cerr << "Processing base_neighbourhood i:" << i << endl;
        // When attacking consider only enemy bases.
        if (base_neighbourhood[i][0]->m_owner == 0)
            continue;

        AttackParameters ap = estimate_attack_feasibility_for_base(base_neighbourhood[i],
                                                                   m_arival_time,
                                                                   max_number_of_attacking_bases);

        if (ap.m_max_arival_time != LOCAL_MAX_INT)
            pq.push(ap);

    }

    //cerr << "Priority queue content:" << endl;
    //while (!pq.empty()) {
    //    cerr << pq.top() << endl;
    //    pq.pop();
    //}

    while(!pq.empty()) {

        int destination_base_id = pq.top().m_destination_base_id;

        if (bases[destination_base_id].m_under_attack == true) {
            pq.pop();
            continue;
        }
        
        if (bases[destination_base_id].m_under_attack == false)
            return AttackParameters(pq.top());
    }

    return AttackParameters();
}


class ExtendedAttackParameters {
public:
    int m_destination_base_id;
    int m_arival_time;
    vector<bool> m_bases_taking_part_in_the_attack;

    ExtendedAttackParameters() { 
        m_destination_base_id = -1; 
        m_arival_time = -1;
        m_bases_taking_part_in_the_attack = vector<bool>();
        }

   ExtendedAttackParameters(int dest, int mat, int n): m_destination_base_id(dest), 
                                                       m_arival_time(mat) {
                                                       m_bases_taking_part_in_the_attack.resize(n, false);
                                                       }

   ExtendedAttackParameters(const ExtendedAttackParameters &eap): m_destination_base_id(eap.m_destination_base_id),
                                                                  m_arival_time(eap.m_arival_time) { 
                                                    //cerr << "Calling AttackParameters copy constructor." << endl;
                                                                  m_bases_taking_part_in_the_attack = eap.m_bases_taking_part_in_the_attack;
                                                                  }
};


void print_eap(ExtendedAttackParameters &eap) {

    cerr << "Printing eap..." << endl;
    cerr << "eap.m_destination_base_id: " << eap.m_destination_base_id << endl; 
    cerr << "eap.m_arival_time: " << eap.m_arival_time << endl;

    cerr << "eap.m_bases_taking_part_in_the_attack: " << endl;
    for(int i = 0; i < eap.m_bases_taking_part_in_the_attack.size(); i++)
        cerr << eap.m_bases_taking_part_in_the_attack[i] << " ";
    cerr << endl;
}



class AbstractWars {
public:
    static int m_step;
    int m_speed;
    int m_no_of_attacking_bases;
    vector<Base> m_bases;
    vector<vector<double>> m_distance_matrix;
    vector<vector<int>> m_arival_time;

    vector<vector<Base*>> m_base_neighbourhood_1;
    vector<vector<Base*>> m_base_neighbourhood_2;
    vector<vector<Base*>> m_base_neighbourhood_3;

    vector<AttackParameters> m_attack_parameters_per_base;

    AbstractWars() {};
    ExtendedAttackParameters single_coordinated_attack(vector<Base*> &base_neighbourhood);
    void prepare_multiple_coordinated_attack(vector<vector<Base*>> &base_neighbourhood);
    void update_bases(vector<int> &bases);
    int init(vector <int> base_locations, int speed);
    vector <int> sendTroops(vector <int> bases, vector <int> troops);
};


ExtendedAttackParameters AbstractWars::single_coordinated_attack(vector<Base*> &base_neighbourhood) {

    int N = base_neighbourhood.size();
    int id = base_neighbourhood[0]->m_id;

    int gathered_troops = 0;
    int max_arival_time = LOCAL_MAX_INT;
    vector<bool> attacking_bases_marker(N, false);

    bool green_light_for_attack = false;

    // Iterate over neighbours of the base.
    for(int i = 1; i < N; i++) {

        if (base_neighbourhood[i]->m_owner != 0)
            continue;

        if (m_bases[i].m_attacking == true)
            continue;

        // How much troops can we send from this base?
        gathered_troops = gathered_troops + (base_neighbourhood[i]->m_size)/2;

        // How long will it take the troops to arrive at the enemies base.
        int attacking_base_id = base_neighbourhood[i]->m_id;
        int troop_arrival_time = m_arival_time[id][attacking_base_id];      

        // What will be the size of the enemies base at troop arival time?
        int future_base_size = estimate_size_in_t_steps(*base_neighbourhood[0], troop_arrival_time);

        //cerr << "troop_arrival_time: " << troop_arrival_time << " gathered_troops: " << gathered_troops << " future_base_size: " << future_base_size << " growth rate: " << base_neighbourhood[0]->m_gr << endl;

        // If there are enough troops attack.
        if (1.2*future_base_size < gathered_troops) {

            green_light_for_attack = true;
            max_arival_time = troop_arrival_time;
            break;
        } else {
            // If there are not enough troops continue gathering troops.
            // Mark the bases that will potentaly attack.
            attacking_bases_marker[i] = true;            
            continue;
        }
    } // for end


    if (green_light_for_attack == true) {

        ExtendedAttackParameters eap(id, max_arival_time, N);
        for(int i = 0; i < N; i++) {
            if (attacking_bases_marker[i] == true) {
                eap.m_bases_taking_part_in_the_attack[i] = true;
                m_bases[i].m_attacking = true;
                m_no_of_attacking_bases++;
            }
        }

        return eap;
    } else {
        return ExtendedAttackParameters();
    }
}


void AbstractWars::prepare_multiple_coordinated_attack(vector<vector<Base*>> &base_neighbourhood) {

    int N = base_neighbourhood.size();

    vector<ExtendedAttackParameters> eap_vec;
    for(int i = 0; i < N; i++) {

        if (m_no_of_attacking_bases == m_bases.size())
            break;

        if (base_neighbourhood[i][0]->m_owner == 0)
            continue;

        cerr << "Enemy cell!" << endl;
        ExtendedAttackParameters eap = single_coordinated_attack(base_neighbourhood[i]);
        print_eap(eap);


    }

}


void AbstractWars::update_bases(vector<int> &bases) {

    int N = bases.size()/2;
    m_no_of_attacking_bases = 0;

    // Calculate the growth rate on the first iteration.
    if (m_step == 1) {

        for(int i = 0; i < N; i++) {

            int gr = bases[2*i + 1] - m_bases[i].m_size;
            m_bases[i].m_owner = bases[2*i];
            m_bases[i].m_size = bases[2*i + 1];
            m_bases[i].m_gr = gr;
            m_bases[i].m_under_attack = false;
            m_bases[i].m_attacking = false;
            m_bases[i].m_attack_time = -1;
        }


        int index1 = 0;
        int index2 = 0;
        int index3 = 0;

        for(int i = 0; i < N; i++) {

            if (m_bases[i].m_gr == 1)
                m_base_neighbourhood_1.push_back(vector<Base*>(N));
            else if (m_bases[i].m_gr == 2)
                m_base_neighbourhood_2.push_back(vector<Base*>(N));
            else
                m_base_neighbourhood_3.push_back(vector<Base*>(N));

            for(int j = 0; j < N; j++) {

                if (m_bases[i].m_gr == 1)
                    m_base_neighbourhood_1[index1][j] = &m_bases[j];
                else if (m_bases[i].m_gr == 2)
                    m_base_neighbourhood_2[index2][j] = &m_bases[j];
                else
                    m_base_neighbourhood_3[index3][j] = &m_bases[j];

            }

            int id0 = m_bases[i].m_id;
            auto cf = [this, &id0](const Base *b1, const Base *b2){ 

                    int id1 = b1->m_id;
                    int id2 = b2->m_id;

                    return m_distance_matrix[id0][id1] < m_distance_matrix[id0][id2]; 

                };


            if (m_bases[i].m_gr == 1) {
                sort(m_base_neighbourhood_1[index1].begin(), m_base_neighbourhood_1[index1].end(), cf);
                index1++;
            } else if (m_bases[i].m_gr == 2) {
                sort(m_base_neighbourhood_2[index2].begin(), m_base_neighbourhood_2[index2].end(), cf);
                index2++;
            } else {
                sort(m_base_neighbourhood_3[index3].begin(), m_base_neighbourhood_3[index3].end(), cf);
                index3++;
            }

        }

    // On the other iterations updata only the owners and sizes.
    // Growth rate is constant.
    } else {

        for(int i = 0; i < N; i++) {

            if ( (m_bases[i].m_under_attack == true) && (m_bases[i].m_attack_time == m_step) ) {
                m_bases[i].m_under_attack = false;
                m_bases[i].m_attack_time = -1;
            }
            m_bases[i].m_attacking = false;
            m_bases[i].m_owner = bases[2*i];
            m_bases[i].m_size = bases[2*i + 1];
        }
    } // else end
}



int AbstractWars::init(vector <int> base_locations, int speed) {

    int N = base_locations.size()/2;

    //cerr << "Number of base locations: " << N << endl;

    m_bases.resize(N);
    m_distance_matrix.resize(N);
    m_arival_time.resize(N);
    m_attack_parameters_per_base.resize(N);
    m_speed = speed;
    m_no_of_attacking_bases = 0;

    for(int i = 0; i < N; i++) {
        m_bases[i].m_id = i;
        m_bases[i].m_x = base_locations[2*i];
        m_bases[i].m_y = base_locations[2*i + 1];
        m_bases[i].m_under_attack = false;
        m_bases[i].m_attacking = false;
        m_distance_matrix[i].resize(N);
        m_arival_time[i].resize(N);
    }

    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {

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

    //cerr << "Print bases: " << endl;
    //print_vector(m_bases);

    //cerr << "Print distance matrix:" << endl;
    //print_matrix(m_distance_matrix);

    return 0;
}


vector<int> AbstractWars::sendTroops(vector <int> bases, vector <int> troops) {

    //cerr << endl << "Step: " << m_step << endl;
    update_bases(bases);
    print_vector(m_bases);

    vector<int> res;
    if (m_step < 1) {
        m_step++;
        return res;
    }


    cerr << "m_base_neighbourhood_3 size: " << m_base_neighbourhood_3.size() << endl;
    int index = 1;
    int id0 = m_base_neighbourhood_3[index][0]->m_id;
    for(int j = 0; j < m_bases.size(); j++) {

        int id1 = m_base_neighbourhood_3[index][j]->m_id;
        cerr << *m_base_neighbourhood_3[index][j] << " d: " << m_distance_matrix[id0][id1] << endl;
    
    }

    cerr << "m_base_neighbourhood_2 size: " << m_base_neighbourhood_2.size() << endl;
    id0 = m_base_neighbourhood_2[index][0]->m_id;
    for(int j = 0; j < m_bases.size(); j++) {

        int id1 = m_base_neighbourhood_2[index][j]->m_id;
        cerr << *m_base_neighbourhood_2[index][j] << " d: " << m_distance_matrix[id0][id1] << endl;
    
    }

    cerr << "m_base_neighbourhood_1 size: " << m_base_neighbourhood_1.size() << endl;
    id0 = m_base_neighbourhood_1[index][0]->m_id;
    for(int j = 0; j < m_bases.size(); j++) {

        int id1 = m_base_neighbourhood_1[index][j]->m_id;
        cerr << *m_base_neighbourhood_1[index][j] << " d: " << m_distance_matrix[id0][id1] << endl;
    
    }


    //ExtendedAttackParameters eap = single_coordinated_attack(m_base_neighbourhood_3[1]);

    
    prepare_multiple_coordinated_attack(m_base_neighbourhood_3);

    /*
    int t_steps = 10;
    int index = 3;
    int fs = estimate_size_in_t_steps(m_bases[index], t_steps);

    //cerr << "current: " << m_bases[index].m_size << " fs: " << fs << endl;

    int orig = 0;
    int dest = 10;

    //cerr << "troop arival time: " << troop_arrival_time(orig, dest, m_distance_matrix, m_speed) << endl;
    //cerr << m_arival_time[orig][dest] << endl;
        
    index = 4;
    for(int j = 0; j < m_bases.size(); j++) {

        int id1 = m_base_neighbourhood[index][j]->m_id;
        //cerr << *m_base_neighbourhood[index][j] << " d: " << m_distance_matrix[index][id1] << endl;
    
    }

    int max_number_of_attacking_bases = int(bases.size()/4);

    AttackParameters ap = estimate_attack_feasibility(m_base_neighbourhood,
                                                      m_arival_time,
                                                      max_number_of_attacking_bases,
                                                      m_bases);


    




    //cerr << "Returned ap:" << endl;
    //cerr << ap << endl;
    if (ap.m_destination_base_id != -1) {

        cerr << "Dispatching troops..." << endl;
        int id = ap.m_destination_base_id;
        int at = ap.m_max_arival_time;
        int gr = ap.m_destination_base_gr;
        cerr << "gr: " << gr << endl;
        
        for(int i = 0; i <= ap.m_boarder_base; i++) {
            if ( m_base_neighbourhood[id][i]->m_owner != 0 )
                continue;

            int from_base = m_base_neighbourhood[id][i]->m_id;
            int to_base = id;

            cerr << "from: " << from_base << " to: " << to_base << endl;

            res.push_back(from_base);
            res.push_back(to_base);

            m_bases[to_base].m_under_attack = true;
            m_bases[to_base].m_attack_time = m_step + at;

        }

    }

    */

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
