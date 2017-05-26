#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <queue>
#include <random>
#include <set>

using namespace std;

const int LOCAL_MAX_INT = numeric_limits<int>::max();
const int MAX_NUMBER_OF_COORDINATED_ATTACKS_PER_STEP_PER_GROWTH_TYPE = 1;
const int MAX_NUMBER_OF_ATTACKING_BASES_IN_COORDINATED_ATTACK = 4;
const int MGR = 3; // max growth rate


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


class BasicAttackParameters {
public:
    int m_origin_base_id;
    int m_destination_base_id;
    int m_arival_time;

    BasicAttackParameters() {
        m_origin_base_id = -1;
        m_destination_base_id = -1;
        m_arival_time = -1;
        }

    BasicAttackParameters(int orig, int dest, int mat, int n): m_origin_base_id(orig),
                                                               m_destination_base_id(dest), 
                                                               m_arival_time(mat) { }

    BasicAttackParameters(const BasicAttackParameters &eap): m_origin_base_id(eap.m_origin_base_id),
                                                             m_destination_base_id(eap.m_destination_base_id),
                                                             m_arival_time(eap.m_arival_time) { }
};


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


void print_eap(const ExtendedAttackParameters &eap) {

    cerr << "Printing eap..." << endl;
    cerr << "eap.m_destination_base_id: " << eap.m_destination_base_id << endl; 
    cerr << "eap.m_arival_time: " << eap.m_arival_time << endl;

    cerr << "eap.m_bases_taking_part_in_the_attack: " << endl;
    for(int i = 0; i < eap.m_bases_taking_part_in_the_attack.size(); i++)
        cerr << eap.m_bases_taking_part_in_the_attack[i] << " ";
    cerr << endl;
}


class PriorityQueueExtendedAttackParametersCompare {
	public:
		inline bool operator ()(const ExtendedAttackParameters &eap1, const ExtendedAttackParameters &eap2) const {
			return eap1.m_arival_time < eap2.m_arival_time;
		}
};


class AbstractWars {
public:

    int B;
    vector<int> baseX, baseY;
    vector<int> others;
    int getRandomBase(int sourceInd);

    default_random_engine m_engine;
    static int m_step;
    int m_speed;
    int m_no_of_attacking_bases;
    vector<double> m_powers;
    vector<Base> m_bases;
    vector<vector<double>> m_distance_matrix;
    vector<vector<int>> m_arival_time;
    int m_controlled_bases;

    vector<vector<Base*>> m_base_neighbourhood;
    vector<vector<Base*>> m_base_neighbourhood_1;
    vector<vector<Base*>> m_base_neighbourhood_2;
    vector<vector<Base*>> m_base_neighbourhood_3;

    vector<vector<vector<Base*>>> m_base_neighbourhood_123;

    vector<AttackParameters> m_attack_parameters_per_base;

    AbstractWars() {};
    bool are_there_any_bases_left(int growth_rate);
    ExtendedAttackParameters single_coordinated_attack(vector<Base*> &base_neighbourhood);
    void prepare_multiple_coordinated_attack(vector<vector<Base*>> &base_neighbourhood,
                                             priority_queue<ExtendedAttackParameters, vector<ExtendedAttackParameters>, PriorityQueueExtendedAttackParametersCompare> &pq);
    bool dispatch_multiple_coordinated_attacks(vector<int> &ret,
                                               priority_queue<ExtendedAttackParameters, vector<ExtendedAttackParameters>, PriorityQueueExtendedAttackParametersCompare> &pq);    
    int get_random_base(int origin);
    void nearest_neighbour_attack(vector<int> &res);
    void attack_random_base(vector<int> &res);
    void refill_zero_base(vector<int> &res);
    void support_bases(vector<int> &res);
    void update_bases(vector<int> &bases, vector<int> &troops);
    int init(vector <int> base_locations, int speed);
    vector <int> sendTroops(vector <int> bases, vector <int> troops);
};


int AbstractWars::getRandomBase(int sourceInd) {
    vector<double> probs(others.size());
    double sp = 0;
    for (int i = 0; i < (int)others.size(); ++i) {
        int ind = others[i];
        probs[i] = 1 / (pow(baseX[sourceInd] - baseX[ind], 2) + pow(baseY[sourceInd] - baseY[ind], 2));
        sp += probs[i];
    }
    
    double r = rand() * 1.0/RAND_MAX * sp;
    double s = 0;
    for (int i = 0; i < (int)others.size(); ++i) {
        s += probs[i];
        if (s >= r)
            return others[i];
    }
    return others[others.size() - 1];
}


bool AbstractWars::are_there_any_bases_left(int growth_rate) {

    for(int i = 0; i < m_base_neighbourhood_123[growth_rate].size(); i++)
        if (m_base_neighbourhood_123[growth_rate][i][0]->m_owner != 0)
            return true;

    return false;

}


ExtendedAttackParameters AbstractWars::single_coordinated_attack(vector<Base*> &base_neighbourhood) {

    int N = base_neighbourhood.size();
    int id = base_neighbourhood[0]->m_id;

    int gathered_troops = 0;
    int max_arival_time = LOCAL_MAX_INT;
    vector<bool> attacking_bases_marker(N, false);

    bool green_light_for_attack = false;

    // Iterate over neighbours of the base.
    int checked_bases = 0;


    double mean = double(m_controlled_bases);
    double std = mean;
    //std::lognormal_distribution<> d(mean, std);
    //std::exponential_distribution<> d(mean);    
    //uniform_int_distribution<> d(log(mean), mean);

    double mean_enemy_distance = 0.0;
    for(int i = 0; i < N; i++) {
        if (base_neighbourhood[i]->m_owner == 0)
            continue;

        int attacking_base_id = base_neighbourhood[i]->m_id;
        mean_enemy_distance += m_distance_matrix[id][attacking_base_id];
    }
    mean_enemy_distance = mean_enemy_distance/(N-1);

    for(int i = 1; i < N; i++) {

        int attacking_base_id = base_neighbourhood[i]->m_id;

        if (base_neighbourhood[i]->m_owner != 0)
            continue;


        if (base_neighbourhood[i]->m_attacking == true)
            continue;

        if(checked_bases > MAX_NUMBER_OF_ATTACKING_BASES_IN_COORDINATED_ATTACK)
            break;

        if(m_distance_matrix[id][attacking_base_id] > 300)
            break;

        checked_bases++;
        
        int origin_id = base_neighbourhood[i]->m_id;

        // How much troops can we send from this base?
        gathered_troops = gathered_troops + (base_neighbourhood[i]->m_size)/2;

        // How long will it take the troops to arrive at the enemies base.
        int troop_arrival_time = m_arival_time[id][attacking_base_id];      

        // What will be the size of the enemies base at troop arival time?
        int future_base_size = estimate_size_in_t_steps(*base_neighbourhood[0], troop_arrival_time);

        // If there are enough troops attack.
        if (1.2*future_base_size < gathered_troops) {
            //cerr << "Firing!!!" << endl;
            green_light_for_attack = true;
            max_arival_time = troop_arrival_time;
            attacking_bases_marker[origin_id] = true;
            break;
        } else {
            // If there are not enough troops continue gathering troops.
            // Mark the bases that will potentaly attack.
            attacking_bases_marker[origin_id] = true;            
            continue;
        }
    } // for end


    if (green_light_for_attack == true) {

        ExtendedAttackParameters eap(id, max_arival_time, N);
        for(int i = 0; i < N; i++) {
            if (attacking_bases_marker[i] == true) {
                eap.m_bases_taking_part_in_the_attack[i] = true;
            }
        }

        return eap;
    } else {
        return ExtendedAttackParameters();
    }
}


void AbstractWars::prepare_multiple_coordinated_attack(vector<vector<Base*>> &base_neighbourhood,
                                                       priority_queue<ExtendedAttackParameters, vector<ExtendedAttackParameters>, PriorityQueueExtendedAttackParametersCompare> &pq) {

    int N = base_neighbourhood.size();
    for(int i = 0; i < N; i++) {

        if (m_no_of_attacking_bases == m_bases.size())
            break;
     
        if (base_neighbourhood[i][0]->m_owner == 0)
            continue;

        if (base_neighbourhood[i][0]->m_under_attack == true)
            continue;

        ExtendedAttackParameters eap = single_coordinated_attack(base_neighbourhood[i]);
        if (eap.m_destination_base_id != -1)
            pq.push(eap);

    }
}


bool AbstractWars::dispatch_multiple_coordinated_attacks(vector<int> &res,
                                                         priority_queue<ExtendedAttackParameters, vector<ExtendedAttackParameters>, PriorityQueueExtendedAttackParametersCompare> &pq) {

    int number_of_coordinated_attacks_per_step_pre_grwoth_type = 0;
    bool return_value = false;

    if (pq.empty() == true)
        return return_value;
    else
        return_value = true;

    while(!pq.empty()) {

        if (number_of_coordinated_attacks_per_step_pre_grwoth_type >= MAX_NUMBER_OF_COORDINATED_ATTACKS_PER_STEP_PER_GROWTH_TYPE)
            break;

        ExtendedAttackParameters eap = pq.top();

        int destination_id = eap.m_destination_base_id;
        m_bases[destination_id].m_under_attack = true;
        m_bases[destination_id].m_attack_time = m_step + eap.m_arival_time;
        for(int i = 0; i < m_bases.size(); i++) {

            if (eap.m_bases_taking_part_in_the_attack[i] == false)
                continue;

            int origin_id = m_bases[i].m_id;
            m_bases[origin_id].m_attacking = true;
            m_no_of_attacking_bases++;

            res.push_back(origin_id);
            res.push_back(destination_id);

        }
        pq.pop();
        number_of_coordinated_attacks_per_step_pre_grwoth_type++;
    }
    return return_value;
}





void AbstractWars::nearest_neighbour_attack(vector<int> &res) {

    int N = m_bases.size();
    for(int i = 0; i < N; i++) {
        
        if (m_base_neighbourhood[i][0]->m_attacking == true)
            continue;

        if (m_base_neighbourhood[i][0]->m_owner != 0)
            continue;

        int origin_id = m_base_neighbourhood[i][0]->m_id;
        int max_arival_time = LOCAL_MAX_INT;
        vector<bool> attacking_bases_marker(N, false);
        bool green_light_for_attack = false;


        // How many troops can we send from this base?
        int gathered_troops = (m_base_neighbourhood[i][0]->m_size)/2;
        for(int j = 1; j < N; j++) {

            if (m_base_neighbourhood[i][j]->m_owner == 0)
                continue;

            // With this if the most distant bases that are full or nearly full
            // will not attack. They will wait and wait.
            // Thus we do not check for m_under_attack.
            // This is an aggresive strategy.
            //if (m_base_neighbourhood[i][j]->m_under_attack == true)
            //    continue;

            // Potential target found!
            green_light_for_attack = true;

            // How long will it take the troops to arrive at the enemies base.
            int attacked_base_id = m_base_neighbourhood[i][j]->m_id;
            int troop_arrival_time = m_arival_time[i][attacked_base_id];      

            // What will be the size of the enemies base at troop arival time?
            int future_base_size = estimate_size_in_t_steps(*m_base_neighbourhood[i][j], troop_arrival_time);

            if (1.2*future_base_size < gathered_troops) {
                res.push_back(origin_id);
                res.push_back(attacked_base_id);

                m_bases[origin_id].m_attacking = true;
                m_bases[attacked_base_id].m_under_attack = true;
                m_bases[attacked_base_id].m_attack_time = m_step + troop_arrival_time;
                break;
            }

            if (green_light_for_attack == true)
                break;
        } // inner for loop

    } // outer for loop


}


void AbstractWars::attack_random_base(vector<int> &res) {

    int N = m_bases.size();
    vector<double> p(N, 0.0);

    bool break_flag = false;

    for(int i = 0; i < N; i++) {

        if (m_bases[i].m_attacking == true)
            continue;

        if (m_bases[i].m_owner != 0)
            continue;
   
        if (m_bases[i].m_size < 900)
            continue;

        double cp = 0.0;
        for(int j = 0; j < N; j++) {

            if (i == j)
                continue;

            if (m_bases[j].m_owner == 0)
                continue;

            double base_p = 1.0/m_distance_matrix[i][j];
            p[j] = cp + base_p;
            cp = cp + base_p;
        }

        uniform_real_distribution<> rd(0.0, 1.0);
        double ap = rd(m_engine)*cp;
        for(int j = 0; j < N; j++) {
            if (p[j] > ap) {
                //cerr << "Random attack!" << endl;
                m_bases[i].m_attacking = true;
                res.push_back(i);
                res.push_back(j);
                break;
            }
        }
    }
}


void AbstractWars::refill_zero_base(vector<int> &res) {

    int N = m_bases.size();
    for(int i = 0; i < N; i++) {

        if (m_bases[i].m_size != 0)
            continue;

        if (m_bases[i].m_under_attack == true)
            continue;

        for(int j = 1; j < N; j++) {

            if (m_base_neighbourhood[i][j]->m_owner != 0)
                continue;

            if (m_base_neighbourhood[i][j]->m_size == 0)
                continue;

            //cerr << "Refilling base!!!" << endl;
            int origin_id = m_base_neighbourhood[i][j]->m_id;
            res.push_back(m_bases[origin_id].m_id);
            res.push_back(m_bases[i].m_id);

            m_bases[origin_id].m_attacking = true;
            m_bases[i].m_under_attack = true;
            m_bases[i].m_attack_time = m_step + m_arival_time[i][origin_id];

            break;

        }
    }
}


void AbstractWars::support_bases(vector<int> &res) {

    int N = m_base_neighbourhood.size();

    for(int i = 0; i < N; i++) {
        
        if (m_base_neighbourhood[i][0]->m_owner != 0)
            continue;

        if (m_base_neighbourhood[i][0]->m_attacking == true)
            continue;

        if (m_base_neighbourhood[i][0]->m_size < 900)
            continue;

        for(int j = 1; j < N; j++) {

            if (m_base_neighbourhood[i][j]->m_attacking == true)
                continue;

            if (m_base_neighbourhood[i][j]->m_owner != 0)
                continue;
            else {

                res.push_back(i);
                res.push_back(m_base_neighbourhood[i][j]->m_id);
                break;
            }
        }
    }
}



void AbstractWars::update_bases(vector<int> &bases, vector<int> &troops) {

    int N = bases.size()/2;

    m_no_of_attacking_bases = 0;
    m_controlled_bases = 0;

    // Calculate the growth rate on the first iteration.
    if (m_step == 0) {

        for(int i = 0; i < N; i++) {

            m_bases[i].m_owner = bases[2*i];

            if(bases[2*i] == 0)
                m_controlled_bases++;

            m_bases[i].m_size = bases[2*i + 1];
            m_bases[i].m_under_attack = false;
            m_bases[i].m_attacking = false;
            m_bases[i].m_attack_time = -1;
        }

   
    } else if (m_step == 1) {

        set<int> player_set;
        for(int i = 0; i < N; i++) {

            int gr = bases[2*i + 1] - m_bases[i].m_size;
            m_bases[i].m_owner = bases[2*i];
            player_set.insert(bases[2*i]);

            if(bases[2*i] == 0)
                m_controlled_bases++;

            m_bases[i].m_size = bases[2*i + 1];
            m_bases[i].m_gr = gr;
            m_bases[i].m_under_attack = false;
            m_bases[i].m_attacking = false;
            m_bases[i].m_attack_time = -1;
        }

        m_powers.resize(player_set.size(), 0);

        int index1 = 0;
        int index2 = 0;
        int index3 = 0;

        vector<int> index_vec(MGR + 1, 0);
        for(int i = 0; i <= MGR; i++)
            m_base_neighbourhood_123.push_back(vector<vector<Base*>>());

        for(int i = 0; i < N; i++) {

            for(int k = 1; k <= MGR; k++)            
                if (m_bases[i].m_gr == k)
                    m_base_neighbourhood_123[k].push_back(vector<Base*>(N));

            if (m_bases[i].m_gr == 1)
                m_base_neighbourhood_1.push_back(vector<Base*>(N));
            else if (m_bases[i].m_gr == 2)
                m_base_neighbourhood_2.push_back(vector<Base*>(N));
            else
                m_base_neighbourhood_3.push_back(vector<Base*>(N));

            for(int j = 0; j < N; j++) {

                for(int k = 1; k <= MGR; k++)
                    if (m_bases[i].m_gr == k)
                        m_base_neighbourhood_123[k][index_vec[k]][j] = &m_bases[j];


                if (m_bases[i].m_gr == 1)
                    m_base_neighbourhood_1[index1][j] = &m_bases[j];
                else if (m_bases[i].m_gr == 2)
                    m_base_neighbourhood_2[index2][j] = &m_bases[j];
                else
                    m_base_neighbourhood_3[index3][j] = &m_bases[j];

                m_base_neighbourhood[i][j] = &m_bases[j];

            }

            int id0 = m_bases[i].m_id;
            auto cf = [this, &id0](const Base *b1, const Base *b2){ 

                    int id1 = b1->m_id;
                    int id2 = b2->m_id;

                    return m_distance_matrix[id0][id1] < m_distance_matrix[id0][id2]; 

                };


            for(int k = 1; k <= MGR; k++) {
                if (m_bases[i].m_gr == k) {
                    sort(m_base_neighbourhood_123[k][index_vec[k]].begin(), m_base_neighbourhood_123[k][index_vec[k]].end(), cf);
                    index_vec[k]++;
                }
            }


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
            sort(m_base_neighbourhood[i].begin(), m_base_neighbourhood[i].end(), cf);
        }

    // On other iterations updata only the owners and sizes.
    // Growth rate is constant.
    } else {

        for(int i = 0; i < m_powers.size(); i++)
            m_powers[i] = 0;

        for(int i = 0; i < N; i++) {

            if ( (m_bases[i].m_under_attack == true) && (m_bases[i].m_attack_time == m_step) ) {
                m_bases[i].m_under_attack = false;
                m_bases[i].m_attack_time = -1;
            }

            if(bases[2*i] == 0)
                m_controlled_bases++;

            m_bases[i].m_attacking = false;
            m_bases[i].m_owner = bases[2*i];
            m_bases[i].m_size = bases[2*i + 1];
            m_powers[ bases[2*i] ] = m_powers[ bases[2*i] ] + bases[2*i + 1];

        }


        int M = troops.size()/4;
        for(int i = 0; i < M; i++)
            m_powers[ troops[4*i] ] = m_powers[ troops[4*i] ] + troops[4*i + 1];

        int tt = 0;
        for(int i = 0; i < m_powers.size(); i++)
            tt = tt + m_powers[i];

        for(int i = 0; i < m_powers.size(); i++)
            m_powers[i] = m_powers[i]/tt;

    } // else end
}


int AbstractWars::init(vector <int> base_locations, int speed) {

    srand(123);
    B = base_locations.size() / 2;
    for (int i = 0; i < B; ++i) {
        baseX.push_back(base_locations[2*i]);
        baseY.push_back(base_locations[2*i+1]);
    }

    random_device rd;
    m_engine.seed(rd());

    int N = base_locations.size()/2;

    m_bases.resize(N);
    m_distance_matrix.resize(N);
    m_arival_time.resize(N);
    m_attack_parameters_per_base.resize(N);
    m_base_neighbourhood.resize(N);
    m_speed = speed;
    m_no_of_attacking_bases = 0;
    m_step = 0;
    m_controlled_bases = 0;

    for(int i = 0; i < N; i++) {
        m_bases[i].m_id = i;
        m_bases[i].m_x = base_locations[2*i];
        m_bases[i].m_y = base_locations[2*i + 1];
        m_bases[i].m_under_attack = false;
        m_bases[i].m_attacking = false;
        m_distance_matrix[i].resize(N);
        m_arival_time[i].resize(N);
        m_base_neighbourhood[i].resize(N);
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

    return 0;
}


vector<int> AbstractWars::sendTroops(vector <int> bases, vector <int> troops) {

    update_bases(bases, troops);
    //print_vector(m_bases);

    vector<int> res;
    if (m_step < 1) {
        m_step++;
        return res;
    }


    //priority_queue<ExtendedAttackParameters, vector<ExtendedAttackParameters>, PriorityQueueExtendedAttackParametersCompare> pq1;
    //prepare_multiple_coordinated_attack(m_base_neighbourhood_123[3], pq1);
    //attack3 = dispatch_multiple_coordinated_attacks(res, pq1);

    //priority_queue<ExtendedAttackParameters, vector<ExtendedAttackParameters>, PriorityQueueExtendedAttackParametersCompare> pq2;
    //prepare_multiple_coordinated_attack(m_base_neighbourhood_123[2], pq2);
    //attack2 = dispatch_multiple_coordinated_attacks(res, pq2);

    //priority_queue<ExtendedAttackParameters, vector<ExtendedAttackParameters>, PriorityQueueExtendedAttackParametersCompare> pq2;
    //prepare_multiple_coordinated_attack(m_base_neighbourhood_123[1], pq2);
    //dispatch_multiple_coordinated_attacks(res, pq2);

    if (m_speed < 7) {
        others.resize(0);
        for (int i = 0; i < B; ++i)
            if (bases[2*i] != 0)
                others.push_back(i);
        if (others.size() == 0) {
            // noone to fight!
            return vector<int>(0);
        }

        for (int i = 0; i < B; ++i) {
            //if (bases[2*i] == 0 && bases[2*i+1] > 1000 * 2 / 3) {
            if (bases[2*i] == 0 && bases[2*i+1] > 975) {
                // send troops to a random base of different ownership
                res.push_back(i);
                res.push_back(getRandomBase(i));
            }
        }
    } else {

        priority_queue<ExtendedAttackParameters, vector<ExtendedAttackParameters>, PriorityQueueExtendedAttackParametersCompare> pq1;
        prepare_multiple_coordinated_attack(m_base_neighbourhood, pq1);
        dispatch_multiple_coordinated_attacks(res, pq1);

    }

    //refill_zero_base(res);

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
