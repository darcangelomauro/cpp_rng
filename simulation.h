#ifndef METROPOLIS_H
#define METROPOLIS_H


class metropolis
{
    public:
        
        // sweep function
        void sweep(workspace w); 
        
        // simulation function
        void simulation(workspace w); 

        // scale tutning function
        void tune_scale(workspace w);

        // print obserables
        void print(workspace w);





    private:

        // INPUT PARAMETERS
        
        // number of sweeps
        int Nsw;

        // number of sweeps between two outputs
        int gap;
        
        // metropolis scale factor
        double scale;

        // END INPUT PARAMETERS

}



#endif
