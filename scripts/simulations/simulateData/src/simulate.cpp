#include <Rcpp.h>
#include <algorithm>
#include <limits>

using namespace Rcpp;

struct Lineage {
    int anc_index;
    int index;
    double start_time;
    double end_time;
    int state;
    std::string status; 
};

// [[Rcpp::export]]
List simulateChromoHiSSEInternal(std::vector<double> lambda,
                                 std::vector<double> mu,
                                 NumericMatrix eta,
                                 List Omega,
                                 double rho,
                                 double time,
                                 int initial_state) {
    
    // initialize the lineages
    typedef std::vector<Lineage> LineageVector;

    // helper stuff
    int num_states = lambda.size();

    // initialize lists of lineages
    LineageVector inactive_lineages;
    LineageVector active_lineages;

    // initialize the two lineages at the root
    Lineage left_lineage;
    left_lineage.anc_index = 0;
    left_lineage.index = 1;
    left_lineage.start_time = time;
    left_lineage.state = initial_state;
    left_lineage.status = "active";

    Lineage right_lineage;
    right_lineage.anc_index = 0;
    right_lineage.index = 2;
    right_lineage.start_time = time;
    right_lineage.state = initial_state;
    right_lineage.status = "active";

    active_lineages.push_back(left_lineage);
    active_lineages.push_back(right_lineage);

    int lineage_counter = 2;

    // simulate
    while ( active_lineages.size() > 0 ) {

        if ( lineage_counter > 50000) {
            Rcout << "maximum number of lineages reached" << std::endl;
            return R_NilValue;
        }

        // get the active lineage to simulate
        Lineage this_lineage = active_lineages.back();
        active_lineages.pop_back();

        // get the current time of the lineage
        double current_time = this_lineage.start_time;

        // Rcout << "simulating lineage " << this_lineage.index << " starting at time " << current_time << std::endl;

        // get the state for the lineage
        int this_state     = this_lineage.state;
        double this_lambda = lambda.at(this_state);
        double this_mu     = mu.at(this_state);
        double this_eta    = -eta(this_state, this_state);

        // simulate this lineage
        bool done = false;
        while (!done) {

            // get current total rate of events
            double total_rate = this_lambda + this_mu + this_eta;

            // draw the next time
            double waiting_time = rexp(1, total_rate)[0];
            current_time -= waiting_time;

            if (current_time < 0.0) {

                // do present sampling
                // Rcout << "rho event" << std::endl;

                // truncate time
                current_time = 0.0;

                // determine if the lineage is sampled
                double v = runif(1, 0, 1)[0];
                if (v < rho) {
                    // we were sampled
                    this_lineage.status = "sampled";
                    // Rcout << "Lineage was sampled at time " << current_time << std::endl;
                } else {
                    // we were not sampled
                    this_lineage.status = "extinct";
                    // Rcout << "Lineage wasn't sampled at time " << current_time << std::endl;
                }

                // update the lineage
                this_lineage.end_time = 0.0;

                // move the lineage to inactive category
                inactive_lineages.push_back(this_lineage);

                // terminate simulation of this lineage
                done = true;

            } else {

                // do event

                // determine the type of event
                double u = runif(1, 0, total_rate)[0];

                if (u < this_lambda) {

                    // speciation event
                    // Rcout << "speciation event" << std::endl;

                    // update this lineage
                    this_lineage.end_time = current_time;
                    this_lineage.status = "split";

                    // move the lineage to inactive category
                    inactive_lineages.push_back(this_lineage);

                    // simulate the daughter states
                    int ancestral_state = this_lineage.state;
                    NumericMatrix this_omega = Omega[ancestral_state];
                    
                    // Rcout << "ancestral state: " << ancestral_state << std::endl;
                    
                    // get left daughter state
                    NumericVector left_daughter_prob = rowSums(this_omega);
                    int left_daughter_state;
                    double v = runif(1, 0, 1)[0];
                    for(int i = 0; i < num_states; ++i) {
                        v -= left_daughter_prob(i);
                        if (v < 0) {
                            left_daughter_state = i;
                            break;
                        }
                    }
                    // Rcout << "simulating left daughter" << std::endl;
                    // Rcout << left_daughter_prob << std::endl;
                    // Rcout << left_daughter_state << std::endl;
                        
                    // get the right daughter state
                    int right_daughter_state;
                    NumericVector right_daughter_prob = this_omega(left_daughter_state, _);
                    v = runif(1, 0, sum(right_daughter_prob))[0];
                    for(int i = 0; i < num_states; ++i) {
                        v -= right_daughter_prob(i);
                        if (v < 0) {
                            right_daughter_state = i;
                            break;
                        }
                    }
                    // Rcout << "simulating right daughter" << std::endl;
                    // Rcout << right_daughter_prob << std::endl;
                    // Rcout << right_daughter_state << std::endl;
                    // stop("STOP");
                    
                    // create new lineages
                    Lineage left_daughter;
                    left_daughter.anc_index = this_lineage.index;
                    left_daughter.index = ++lineage_counter;
                    left_daughter.start_time = current_time;
                    left_daughter.state = left_daughter_state;
                    left_daughter.status = "active";

                    Lineage right_daughter;
                    right_daughter.anc_index = this_lineage.index;
                    right_daughter.index = ++lineage_counter;
                    right_daughter.start_time = current_time;
                    right_daughter.state = right_daughter_state;
                    right_daughter.status = "active";

                    // add new lineages to queue
                    active_lineages.push_back(left_daughter);
                    active_lineages.push_back(right_daughter);

                    // terminate simulation of this lineage
                    done = true;

                    // Rcout << "Lineage split at time " << current_time << std::endl;

                } else if (u < this_lambda + this_mu) {

                    // extinction event
                    // Rcout << "extinction event" << std::endl;

                    // update this lineage
                    this_lineage.end_time = current_time;
                    this_lineage.status = "extinct";

                    // move the lineage to inactive category
                    inactive_lineages.push_back(this_lineage);

                    // terminate simulation of this lineage
                    done = true;

                    // Rcout << "Lineage went extinct at time " << current_time << std::endl;

                } else {

                    // state-change event
                    // Rcout << "mutation event" << std::endl;

                    // get the vector of rates away from this lineage
                    NumericVector these_rates = eta(this_state, _);

                    // choose the new state
                    double v = runif(1, 0, this_eta)[0];

                    int new_state = this_state;
                    for(int i = 0; i < num_states; ++i) {
                        if (i != this_state) {
                            v -= these_rates(i);
                            if (v < 0) {
                                new_state = i;
                                break;
                            }
                        }
                    }

                    // update the lineage
                    this_lineage.state = new_state;

                    // update the rates
                    int old_state = this_state;
                    this_state    = new_state;
                    this_lambda   = lambda.at(this_state);
                    this_mu       = mu.at(this_state);
                    this_eta      = -eta(this_state, this_state);

                    // Rcout << "Lineage mutated at time " << current_time << " -- " << old_state << " -- " << new_state << std::endl;

                }

            }

        } // done simulating this lineage

        // Rcout << "done simulating lineage" << std::endl;

        // stop("done simulating lineage");

    }

    // done simulating lineages

    // create the return value

    // reformat for ape
    int num_inactive_lineages = inactive_lineages.size();
    std::vector<int> anc(num_inactive_lineages);
    std::vector<int> desc(num_inactive_lineages);
    std::vector<double> start_time(num_inactive_lineages);
    std::vector<double> end_time(num_inactive_lineages);
    std::vector<int> state(num_inactive_lineages);
    std::vector<std::string> status(num_inactive_lineages);
    for(int i = 0; i < num_inactive_lineages; ++i) {

        // get the lineage
        Lineage& this_lineage = inactive_lineages.at(i);

        // store info
        anc.at(i) = this_lineage.anc_index;
        desc.at(i) = this_lineage.index;
        start_time.at(i) = this_lineage.start_time;
        end_time.at(i) = this_lineage.end_time;
        state[i] = this_lineage.state;
        status.at(i) = this_lineage.status;

    }

    // populate list
    List list = List::create(
        Named("anc") = anc,
        Named("desc") = desc,
        Named("start_time") = start_time,
        Named("end_time") = end_time,
        Named("state") = state,
        Named("status") = status
    );

    // return list
    return list;

    // List list;
    // return list;
    
}