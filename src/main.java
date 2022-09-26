import java.io.*;
import java.util.*;
import java.lang.*;
import java.math.*;

public class main {
    private static int N = 100; //population size; MUST BE DIVISIBLE BY 4 (mating reasons)
    private static double p = 0.05; //mutation rate (per node)
    private static int G = 1000; //number of generations
    private static int S = 100; //number of species; i.e. number of runs of G generations
    private static int rD = 100; //reference distance, for easier fitness scaling. MAKE SURE it is impossible to have a sequence with under this distance
    private static int analysis = 1; //-2 if test, -1 if analyzing sequences in generations, 0 if analyzing generations, 1 if analyzing species

    private static int nodes; //number of destinations (cities)

    private static int[][] a; //matrix of distances
    private static int[][] seqs; //sequences (population)

    private static double bestDistance = Integer.MAX_VALUE;
    private static int[] bestSequence;
    private static int bestGen;
    private static double bestGenDistance = Integer.MAX_VALUE;
    private static int[] bestGenSequence;

    //creates initial sequences (zeroth generation)
    public static int[] newSeq(){
        ArrayList<Integer> thing = new ArrayList<>();
        int[] out = new int[nodes+1];
        for(int i=1; i<nodes; i++){
            thing.add(i);
        }
        for(int i=nodes-1; i>0; i--){
            int r = (int)(Math.random()*i);
            out[nodes-i]=thing.get(r);
            thing.remove(r);
        }
        out[0] = 0;
        out[nodes] = 0;

        return out;
    }

    public static double fitness1(double distance){
        distance = distance/rD;
        return (1.0/(distance));
    }

    public static double fitness2(double distance){
        distance = distance/rD;
        return (1.0/Math.pow(distance,2.0));
    }

    public static double fitness3(double distance){
        distance = distance/rD;
        return (1.0/Math.pow(distance,3.0));
    }

    public static double fitness5(double distance){
        distance = distance/rD;
        return (1.0/Math.pow(distance,5.0));
    }

    public static double fitness10(double distance){
        distance = distance/rD;
        return (Math.pow(10.0,10.0)/Math.pow(distance,10.0));
    }

    public static double fitnessLog100(double distance){
        distance = distance/rD;
        return (Math.pow(10.0,10.0)/Math.pow(Math.log(distance),10.0));
    }

    //new generation process, how to select parents based on fitness
    public static int[][] nextGeneration(int[][] gen1, int gen){
        //initialization
        int[][] gen2 = new int[N][nodes+1];
        int[] distances = new int[N]; //distances of each seequence

        bestGenDistance = Integer.MAX_VALUE;

        //distance assignment to each sequence
        for(int i=0; i<N; i++){
            int distance = 0;
            int[] seq = gen1[i];
            for(int j=0; j<nodes; j++){
                distance+=a[seq[j]][seq[j+1]];
            }
            distances[i]=distance;
            if(distance<bestDistance){
                bestDistance = distance;
                bestSequence = seq;
                bestGen = gen;
            }
            if(distance<bestGenDistance){
                bestGenDistance = distance;
                bestGenSequence = seq;
            }
        }

        int[] parents = parentsRoulette(distances);

        //mating process
        for(int i=0; i<N/4; i++){
            int p1 = parents[2*i];
            int p2 = parents[2*i+1];

            //create 4 offspring
            for(int j=0; j<4; j++){
                gen2[4*i+j]=mutateFrameshift(partialMapCross(gen1[p1],gen1[p2]));
            }
        }
        return gen2;
    }

    public static int[] parentsRoulette(int[] distances){ //roulette selection: probability of selecting is fitness/total fitness
        double[] fitness = new double[N]; //the fitness of each sequence. Fitness = 1/distance; lesser distance, better fitness
        double F = 0; //Total fitness

        for(int i=0; i<N; i++){
            fitness[i]=fitness3(distances[i]);
            F+=fitness[i];
        }

        //creation of list from 1 to N to keep track of unselected sequences
        int[] thing = new int[N];
        for(int i=0; i<N; i++){
            thing[i]=i;
        }

        //selection of N/2 parents for next generation
        int[] parents = new int[N/2]; //index of parents
        for(int i=0; i<N/2; i++){
            double r = Math.random()*F;
            int j=0;
            while(r>0){
                r-=fitness[thing[j]];
                j++;
            }
            parents[i]=thing[j-1];
            F-=fitness[thing[j-1]];
        }

        //        for(int i=0; i<N; i++){
//            System.out.print("Fitness "+fitness[i]+" of: ");
//            for(int j=0; j<nodes+1; j++){
//                System.out.print(gen1[i][j]+" ");
//            }
//            System.out.println();
//        }
//
//        for(int i:parents){
//            System.out.print("Selected: Fitness "+fitness[i]+" of: ");
//            for(int j=0; j<nodes+1; j++){
//                System.out.print(gen1[i][j]+" ");
//            }
//            System.out.println();
//        }

        return parents;
    }

    //crossing over (mating process)
    public static int[] partialMapCross(int[] p1, int[] p2){
        //Partial Mapping Crossover: chooses a random subsequence of a random parent, offspring receives subsequence in the same position, and the rest is filled by other parent with mapping
        //Choose 2 random numbers 0 to nodes-1, that's the subsequence indices
        int[] offspring = new int[nodes+1];
        offspring[0] = 0;
        offspring[nodes] = 0;

        //choosing parent that gives subsequence
        double r = Math.random();
        if(r<0.5){
            p1 = p1;
            p2 = p2;
        }
        else{
            int[] copy = p1;
            p1 = p2;
            p2 = copy;
        }

        int i1 = (int) Math.random()*(nodes-1)+1;
        int i2 = (int) Math.random()*(nodes-1)+1;
        //assigning offspring from one parent i1-i2 otherwise from other parent with mapping
        if(i2>i1){
            i1 = i1;
            i2 = i2;
        }
        else{
            int copy = i1;
            i1 = i2;
            i2 = copy;
        }

        for(int i=1; i<nodes; i++){
            if(i>=i1&&i<=i2){
                offspring[i]=p1[i];
            }
            else{
                boolean thing = false;
                int index = i;
                while(!thing){
                    thing = true;
                    for(int j=i1; j<=i2; j++){
                        if(p1[j]==p2[index]){
                            thing = false;
                            index = j;
                            break;
                        }
                    }
                }
                offspring[i]=p2[index];
            }
        }
        return offspring;
    }

    //mutations
    public static int[] mutateConsecutive(int[] seq){
        //possible to switch indices 1&2,2&3,3&4,...nodes-2&nodes-1 with equal probability
        for(int i=1; i<nodes-1; i++){
            double r = Math.random();
            if(r<p){
                int copy = seq[i];
                seq[i]=seq[i+1];
                seq[i+1]=copy;
            }
        }
        return seq;
    }

    public static int[] mutateRandom(int[] seq){
        //For each node, if probability p, then another node is selected and they are switched
        for(int i=1; i<nodes-1; i++){
            double r = Math.random();
            if(r<p){
                int other = (int)(Math.random()*(nodes-3)+1);
                if(other>=i){
                    other++;
                }
                int copy = seq[i];
                seq[i]=seq[other];
                seq[other]=copy;
            }
        }
        return seq;
    }

    public static int[] mutateFrameshift(int[] seq){
        //For each node, if probability p, then another position is selected, and that node is moved to that position, and all nodes between are moved one up or down to fill in the gap
        for(int i=1; i<nodes-1; i++){
            double r = Math.random();
            if(r<p){
                int other = (int)(Math.random()*(nodes-3)+1);
                //case where we're moving node up, thus frameshifting down
                if(other>=i){
                    other++;
                    int copy = seq[i];
                    for(int j=i; j<other; j++){
                        seq[j]=seq[j+1];
                    }
                    seq[other]=copy;
                }
                //case where we're moving node down, thus frameshifting up
                else{
                    int copy = seq[i];
                    for(int j=i; j>other; j--){
                        seq[j]=seq[j-1];
                    }
                    seq[other]=copy;
                }
            }
        }
        return seq;
    }

    public static int[] mutateInversion(int[] seq){
        //For each node, if probability p, then another node is selected, the sequence between and containing the nodes are inversed
        for(int i=1; i<nodes-1; i++){
            double r = Math.random();
            if(r<p){
                int first = i;
                int other = (int)(Math.random()*(nodes-3)+1);
                if(other>i){
                    first = first;
                    other = other;
                }
                else{
                    int copy = first;
                    first = other;
                    other = copy;
                }
            }
        }
        return seq;
    }

    public static void main(String[] args) throws IOException {
        //initialization
        BufferedReader file = new BufferedReader(new FileReader("in1"));
        String line = file.readLine();
        nodes = Integer.parseInt(line);

        //setting distances
        a = new int[nodes][nodes];
        for(int i=0; i<nodes; i++){
            line = file.readLine();
            String[] thing = line.split(" ");
            for(int j=0; j<nodes; j++){
                a[i][j] = Integer.parseInt(thing[j]);
            }
        }

//        System.out.println("Generation 0: ");
////        for(int i=0; i<N; i++){
////            for(int j=0; j<nodes+1; j++){
////                System.out.print(seqs[i][j]+" ");
////            }
////            System.out.println();
////        }

        if(analysis==-2){ //mutation testing
            int thing[] = new int[nodes+1];
            thing[0] = 0;
            thing[nodes] = 0;
            for(int i=1; i<nodes; i++){
                thing[i]=i;
            }
            for(int i=1; i<=1000; i++){
                thing = mutateFrameshift(thing);
                System.out.print("After "+i+" mutations: ");
                for(int j:thing){
                    System.out.print(j+" ");
                }
                System.out.println();
            }
        }

        else if(analysis==0||analysis==-1) {

            //creating zeroth generation
            seqs = new int[N][nodes+1];
            for(int i=0; i<N; i++){
                seqs[i]=newSeq();
            }

            //running generations
            for (int k = 0; k < G; k++) {
                seqs = nextGeneration(seqs, k);
                System.out.println("Best Distance of Generation " + (k) + ": " + bestGenDistance);
                System.out.print("Best Sequence of Generation " + (k) + ": ");
                for (int i : bestGenSequence) {
                    System.out.print(i + " ");
                }
                System.out.println();
                if(analysis==-1) {
                    System.out.println("Generation " + (k + 1) + ": ");
                    for (int i = 0; i < N; i++) {
                        for (int j = 0; j < nodes + 1; j++) {
                            System.out.print(seqs[i][j] + " ");
                        }
                        System.out.println();
                    }
                }
            }

            System.out.println("Best Overall Distance in Generation " + (bestGen) + " with Distance " + bestDistance);
            System.out.print("Corresponding Best Sequence: ");
            for (int i : bestSequence) {
                System.out.print(i + " ");
            }
        }

        else if(analysis==1) {

            double bestDistanceTotal = 0;
            double bestGenTotal = 0;
            double bestSpeciesDistance = Integer.MAX_VALUE;
            int[] bestSpeciesSequence = new int[nodes+1];

            for (int l = 0; l < S; l++) {

                //creating zeroth generation
                seqs = new int[N][nodes+1];
                for(int i=0; i<N; i++){
                    seqs[i]=newSeq();
                }

                for (int k = 0; k < G; k++) {
                    seqs = nextGeneration(seqs, k);
                }

                bestDistanceTotal+=bestDistance;
                bestGenTotal+=bestGen;

                if(bestDistance<bestSpeciesDistance){
                    bestSpeciesDistance = bestDistance;
                    bestSpeciesSequence = bestSequence;
                }

                bestDistance = Integer.MAX_VALUE;
                bestGenDistance = Integer.MAX_VALUE;
            }

            System.out.println("In "+S+" species, Average Best Distance in "+G+" Generations and "+N+" Sequences in each Generation is "+bestDistanceTotal/S);
            System.out.println("And the Average Generation the Best Distance was achieved in is "+bestGenTotal/S);
            System.out.println("The Overall Best Distance found was "+bestSpeciesDistance);
            System.out.print("Corresponding Best Sequence: ");
            for (int i : bestSpeciesSequence) {
                System.out.print(i + " ");
            }
            System.out.println();

            System.out.println("ABD = "+(int)(bestDistanceTotal/S)+", AG = "+(int)(bestGenTotal/S)+", OBD = "+(int)bestSpeciesDistance);
        }
    }
}
