// geneDrop_RD.c
// C script to perform gene dropping
// compile: gcc geneDrop_final.c -lgsl
// run: a.out [dropfile] [filter ungenotyped individuals? 0 or 1] [data or simulation? D or R] [output type L or S] [XY, ZW, or autosome data? A, Z, or X]

// include packages
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

// macros for easier code reading
#define ALLELE(x,y,z) allele[o3d(x, y, z)]
#define ALLELE_ONE_RUN(x,y) alleleOneRun[o2d(x,y)]
#define ALLELE_ONE_RUN_DATA(x,y) alleleOneRunData[o2d(x,y)]

// individual structure
struct ind {
  int *a1; // you get this gene from mom
  int *a2; // you get this gene from dad
  struct ind *mom;
  struct ind *dad;
  int i;
  int *cohort;
  int geno;
  int cohortNumber;
  int sex;
};

// initialize a bunch of variables
struct ind *tree;
int treeLength; // number of individuals in tree
int *a;
int *cohortIndicator;
int *allele;
int *alleleOneRun;
int *alleleOneRunData;
int *print_filter;

int nIndividuals, nAlleles, nCohorts;

int o3d(int x, int y, int z) { 
    return (z * nCohorts * nAlleles) + (y * nCohorts) + x; 
}

int o2d(int x, int y) {
    return (y * nCohorts) + x;
}

void fill_tree(char *line, int *print_filter, int print_filter_switch) {
    int i, ia1, ia2, imom, idad, icN, isex, j;
    char *ic;

    i = atoi(strtok(line,","));
    tree[i].i = i;

    ia1 = atoi(strtok(NULL,","));
    tree[i].a1 = &a[ia1];

    ia2 = atoi(strtok(NULL,","));
    tree[i].a2 = &a[ia2];

    idad = atoi(strtok(NULL,","));
    if ( idad == -1 ) { tree[i].dad = NULL; } else { tree[i].dad = &tree[idad]; }
    
    imom = atoi(strtok(NULL,","));
    if ( imom == -1 ) { tree[i].mom = NULL; } else { tree[i].mom = &tree[imom]; }

    isex = atoi(strtok(NULL,","));
    tree[i].sex = isex;

    icN = atoi(strtok(NULL,","));
    tree[i].cohortNumber = icN;

    if (icN > 0) {
        j = 0;
        tree[i].cohort = malloc(sizeof(int) * icN);
        ic = strtok(NULL, ",");
        while (ic != NULL) {
            tree[i].cohort[j] = atoi(ic);
            ic = strtok(NULL, ",");
            j = j + 1;
        }
    }

    // setting up print_filter for this individual.
    // if print_filter == 1, print this ind's genotype in output
    // only counting genotyped, non-founder individuals
    if (print_filter_switch == 1) {
        if ( idad != -1 && imom != -1 && ia1 != 0 && ia2 != 0)
        {
            print_filter[i] = 1;
        } else {
            print_filter[i] = 0;
        }
    } else {
        print_filter[i] = 1;
    }
}

// drop_genes_faster_auto function to perform gene dropping for autosomal and pseudoautosomal loci
void drop_genes_faster_auto (gsl_rng *r) {
    // initialize variables
    int i, rnd;
    // for every individual in the tree:
    for (i = 0; i < treeLength; i++) {
    // if mom has a genotype and dad has a genotype:
    if (tree[i].mom != NULL && tree[i].dad != NULL) {
        // roll a 4-sided die
        rnd = gsl_rng_uniform_int(r, 4);
    // if a 0 was rolled:
    if (rnd == 0) {
        // take allele 1 from mom
        tree[i].a1 = tree[i].mom->a1;
        // take allele 1 from dad
        tree[i].a2 = tree[i].dad->a1;
    }
    // if a 1 was rolled:
    if (rnd == 1) {
        // take allele 2 from mom
        tree[i].a1 = tree[i].mom->a2;
        // take allele 1 from dad
        tree[i].a2 = tree[i].dad->a1;
    }
    // if a 2 was rolled:
    if (rnd == 2) {
        // take allele 1 from mom
        tree[i].a1 = tree[i].mom->a1;
        // take allele 2 from dad
        tree[i].a2 = tree[i].dad->a2;
    }
    // if a 3 was rolled:
    if (rnd == 3) {
        // take allele 2 from mom
        tree[i].a1 = tree[i].mom->a2;
        // take allele 2 from dad
        tree[i].a2 = tree[i].dad->a2;
      }
    }
  }
}

// drop_genes_faster_ZW function to perform gene dropping for sex chromosome loci (non-pseudoautosomal) in a ZW system
// this has different transmission rules for males and females
void drop_genes_faster_ZW (gsl_rng *r) {
    // initialize variables
    int i, rnd;
    // for every individual in the tree:
    for (i = 0; i < treeLength; i++) {
    // if mom has a genotype and dad has a genotype:
    if (tree[i].mom != NULL && tree[i].dad != NULL) {
        // flip a coin (or roll a 2-sided die)
        rnd = gsl_rng_uniform_int(r, 2);
    // if a 0 was rolled:
    if (rnd == 0) {
        // Inheritance from mom depends on whether this is a son or a daughter
        // Sex is encoded numerically so 1 = male, 2 = female
        // if this individual is male:
        if (tree[i].sex == 1) {
            // take allele 1 from dad
            tree[i].a2 = tree[i].dad->a1;
            // take allele 2 from mom (a1 is a placeholder in females)
            tree[i].a1 = tree[i].mom->a2;
        }
        // if this individual is female:
        if (tree[i].sex == 2) {
            // take allele 1 from dad
            tree[i].a2 = tree[i].dad->a1;
            // take allele 1 from dad again as a placeholder
            tree[i].a1 = tree[i].dad->a1;
        }
    }
    // if a 1 was rolled:
    if (rnd == 1) {
        // if this individual is male:
        if (tree[i].sex == 1) {
            // take allele 2 from dad
            tree[i].a2 = tree[i].dad->a2;
            // take allele 2 from mom (a1 is a placeholder in females)
            tree[i].a1 = tree[i].mom->a2;
        }
        // if this individual is female:
        if (tree[i].sex == 2) {
            // take allele 2 from dad
            tree[i].a2 = tree[i].dad->a2;
            // take allele 2 from dad again as a placeholder
            tree[i].a1 = tree[i].dad->a2;
        }
      }
    }
  }
}

// drop_genes_faster_XY function to perform gene dropping for sex chromosome loci (non-pseudoautosomal) in an XY system
// this has different transmission rules for males and females
void drop_genes_faster_XY (gsl_rng *r) {
    // initialize variables
    int i, rnd;
    // for every individual in the tree:
    for (i = 0; i < treeLength; i++) {
    // if mom has a genotype and dad has a genotype:
    if (tree[i].mom != NULL && tree[i].dad != NULL) {
        // flip a coin (or roll a 2-sided die)
        rnd = gsl_rng_uniform_int(r, 2);
    // if a 0 was rolled:
    if (rnd == 0) {
        // Inheritance from dad depends on whether this is a son or a daughter
        // Sex is encoded numerically so 1 = male, 2 = female
        // if this individual is female:
        if (tree[i].sex == 2) {
            // take allele 1 from dad (a2 is a placeholder in males)
            tree[i].a2 = tree[i].dad->a1;
            // take allele 1 from mom 
            tree[i].a1 = tree[i].mom->a1;
        }
        // if this individual is male:
        if (tree[i].sex == 1) {
            // take allele 1 from mom
            tree[i].a1 = tree[i].mom->a1;
            // take allele 1 from mom again as a placeholder
            tree[i].a2 = tree[i].mom->a1;
        }
    }
    // if a 1 was rolled:
    if (rnd == 1) {
        // if this individual is female:
        if (tree[i].sex == 2) {
            // take allele 1 from dad (a2 is a placeholder in males)
            tree[i].a2 = tree[i].dad->a1;
            // take allele 2 from mom 
            tree[i].a1 = tree[i].mom->a2;
        }
        // if this individual is male:
        if (tree[i].sex == 1) {
            // take allele 2 from mom 
            tree[i].a1 = tree[i].mom->a2;
            // take allele 2 from mom again as a placeholder
            tree[i].a2 = tree[i].mom->a2;
        }
      }
    }
  }
}

void init_freq() {
    int i, j, k;
    for (i = 0; i < nCohorts; i++) {
        for (j = 0; j < nAlleles; j++) {
            for (k = 0; k < treeLength + 1; k++) {
                ALLELE(i,j,k) = 0;
            }
        }
    }
}

void init_freq_one_run() {
    int i, j;
    for (i = 0; i < nCohorts; i++) {
        for (j = 0; j < nAlleles; j++) {
            ALLELE_ONE_RUN(i,j) = 0;
        }
    }
}

void print_freq_one_run_v2(int run, int allele) {
    int i, j;
    printf("%d ", run);
    for (i = 0; i < nCohorts; i++) {
        if (cohortIndicator[i] == 1) {
            for (j = 0; j < nAlleles; j++) {
                if (j == allele) {
                    printf("%d ", ALLELE_ONE_RUN(i,j));
                }
            }
        }
        else {
            printf("NA ");
        }
    }
    printf("\n");
}
void print_freq_one_run(int run) {
    int i, j, s;
    for (i = 0; i < nCohorts; i++) {
        if (cohortIndicator[i] == 1) {
            s = 0;
            j = 0;
            while (s == 0 && j < nAlleles) {
                s = s + ALLELE_ONE_RUN(i,j);
                j += 1;
            }
            if (s > 0) {
                for (j = 1; j < nAlleles; j++) { /* do not print allele 0 (missing data) */
                    printf("%d %d %d %d\n", run, i, j, ALLELE_ONE_RUN(i,j));
                }
            }
        }
    }
}

// calc_freq_one_run_auto function to calculate the counts of each allele in each cohort for autosomal and pseudoautosomal loci
void calc_freq_one_run_auto(int startingCohort, int * print_filter) {
    // initialize variables
    int i, j, k, l, t;

    // for every individual in the tree:
    for (i = 0; i < treeLength; i++)
    {
        // if we are including this individual:
        if (print_filter[i] == 1)
        {
            // for every allele at this locus:
            for (l = 0; l < nAlleles; l++)
            {
                // number of copies of focal allele = 0
                k = 0;
                // if we have data from the year this individual was born:
                if (tree[i].cohortNumber > 0)
                {
                    // if this individual's 1st allele is the focal allele and the 2nd isn't:
                    if ((*(tree[i].a1) == l) && (*(tree[i].a2) != l)) {
                        // number of copies of this allele = 1
                        k = 1;
                    }
                    // if this individual's 2nd allele is the focal allele and the 1st isn't:
                    else if ((*(tree[i].a1) != l) && (*(tree[i].a2) == l)) {
                        // number of copies of this allele = 1
                        k = 1;
                    }
                    // if both of this individual's alleles are the focal allele:
                    else if ((*(tree[i].a1) == l) && (*(tree[i].a2) == l)) {
                        // number of copies of this allele = 2
                        k = 2;
                    }
                    // for every cohort:
                    for (j = 0; j < tree[i].cohortNumber; j++)
                    {
                        // calculate how many cohorts after the first this is
                        t = tree[i].cohort[j] - startingCohort;
                        // add the number of copies of this allele to the entry for this cohort and this allele in the alleleOneRun table
                        ALLELE_ONE_RUN(t,l) += k;
//                      printf("copy_num: %d ind: %d allele: %d a1: %d a2: %d cohort: %d allele_one_run: %d %d %d-\n", k, i, l, *tree[i].a1, *tree[i].a2, tree[i].cohort[j], ALLELE_ONE_RUN(t,l),t,l);
                    }
                }
            }
        }
    }
}

// calc_freq_one_run_ZW function to calculate the counts of each allele in each cohort for sex chromosome loci in a ZW system
void calc_freq_one_run_ZW(int startingCohort, int * print_filter) {
    // initialize variables
    int i, j, k, l, t;

    // for every individual in the tree:
    for (i = 0; i < treeLength; i++)
    {
        // if we are including this individual:
        if (print_filter[i] == 1)
        {
            // for every allele at this locus:
            for (l = 0; l < nAlleles; l++)
            {
                // number of copies of focal allele = 0
                k = 0;
                // if we have data from the year this individual was born:
                if (tree[i].cohortNumber > 0)
                {
                    // if this individual is male:
                    if (tree[i].sex == 1)
                    {
                        // if this individual's 1st allele is the focal allele and the 2nd isn't:
                        if ((*(tree[i].a1) == l) && (*(tree[i].a2) != l)) {
                            // number of copies of this allele = 1
                            k = 1;
                        }
                        // if this individual's 2nd allele is the focal allele and the 1st isn't:
                        else if ((*(tree[i].a1) != l) && (*(tree[i].a2) == l)) {
                            // number of copies of this allele = 1
                            k = 1;
                        }
                        // if both of this individual's alleles are the focal allele:
                        else if ((*(tree[i].a1) == l) && (*(tree[i].a2) == l)) {
                            // number of copies of this allele = 2
                            k = 2;
                        }
                    }
                    // if this individual is female:
                    else if (tree[i].sex == 2)
                    {
                        // if this individual's 2nd allele is the focal allele:
                        if (*(tree[i].a2) == l) {
                            // number of copies of this allele = 1
                            k = 1;
                        }
                    }
                    // for every cohort:
                    for (j = 0; j < tree[i].cohortNumber; j++)
                    {
                        // calculate how many cohorts after the first this is
                        t = tree[i].cohort[j] - startingCohort;
                        // add the number of copies of this allele to the entry for this cohort and this allele in the alleleOneRun table
                        ALLELE_ONE_RUN(t,l) += k;
                    }
                }
            }
        }
    }
}

// calc_freq_one_run_XY function to calculate the counts of each allele in each cohort for sex chromosome loci in an XY system
void calc_freq_one_run_XY(int startingCohort, int * print_filter) {
    // initialize variables
    int i, j, k, l, t;

    // for every individual in the tree:
    for (i = 0; i < treeLength; i++)
    {
        // if we are including this individual:
        if (print_filter[i] == 1)
        {
            // for every allele at this locus:
            for (l = 0; l < nAlleles; l++)
            {
                // number of copies of focal allele = 0
                k = 0;
                // if we have data from the year this individual was born:
                if (tree[i].cohortNumber > 0)
                {
                    // if this individual is female:
                    if (tree[i].sex == 2)
                    {
                        // if this individual's 1st allele is the focal allele and the 2nd isn't:
                        if ((*(tree[i].a1) == l) && (*(tree[i].a2) != l)) {
                            // number of copies of this allele = 1
                            k = 1;
                        }
                        // if this individual's 2nd allele is the focal allele and the 1st isn't:
                        else if ((*(tree[i].a1) != l) && (*(tree[i].a2) == l)) {
                            // number of copies of this allele = 1
                            k = 1;
                        }
                        // if both of this individual's alleles are the focal allele:
                        else if ((*(tree[i].a1) == l) && (*(tree[i].a2) == l)) {
                            // number of copies of this allele = 2
                            k = 2;
                        }
                    }
                    // if this individual is male:
                    else if (tree[i].sex == 1)
                    {
                        // if this individual's 1st allele is the focal allele:
                        if (*(tree[i].a1) == l) {
                            // number of copies of this allele = 1
                            k = 1;
                        }
                    }
                    // for every cohort:
                    for (j = 0; j < tree[i].cohortNumber; j++)
                    {
                        // calculate how many cohorts after the first this is
                        t = tree[i].cohort[j] - startingCohort;
                        // add the number of copies of this allele to the entry for this cohort and this allele in the alleleOneRun table
                        ALLELE_ONE_RUN(t,l) += k;
                    }
                }
            }
        }
    }
}

// calc_freq_one_run_data_auto function to calculate the number of alleles for each cohort for autosomal and pseudoautosomal loci
// This is needed to calculate allele frequency for data
void calc_freq_one_run_data_auto(int startingCohort, int * print_filter) {
    // initialize variables
    int i, j, k, l, t;

    // for every individual in the tree:
    for (i = 0; i < treeLength; i++)
    {
        // if we're including this individual
        if (print_filter[i] == 1)
        {
            // for every allele at this locus
            for (l = 0; l < nAlleles; l++)
            {
                // number of copies of this allele = 0
                k = 0;
                // if we have data from the year this individual was born:
                if (tree[i].cohortNumber > 0)
                {
                    // for every cohort:
                    for (j = 0; j < tree[i].cohortNumber; j++)
                    {
                        // calculate how many cohorts after the first this is
                        t = tree[i].cohort[j] - startingCohort;
                        // add 2 to the entry for this cohort and this allele in the alleleOneRun table (because this individual has 2 alleles total)
                        ALLELE_ONE_RUN_DATA(t,l) += 2;
                    }
                }
            }
        }
    }
}

// calc_freq_one_run_data_ZW function to calculate the number of alleles for each cohort for sex chromosome loci in a ZW system
// This is needed to calculate allele frequency for data
void calc_freq_one_run_data_ZW(int startingCohort, int * print_filter) {
    // initialize variables
    int i, j, k, l, t;

    // for every individual in the tree:
    for (i = 0; i < treeLength; i++)
    {
        // if we're including this individual
        if (print_filter[i] == 1)
        {
            // for every allele at this locus
            for (l = 0; l < nAlleles; l++)
            {
                // number of copies of this allele = 0
                k = 0;
                // if we have data from the year this individual was born:
                if (tree[i].cohortNumber > 0)
                {
                    // if this individual is male:
                    if (tree[i].sex == 1)
                    {
                        // for every cohort:
                        for (j = 0; j < tree[i].cohortNumber; j++)
                        {
                            // calculate how many cohorts after the first this is
                            t = tree[i].cohort[j] - startingCohort;
                            // add 2 to the entry for this cohort and this allele in the alleleOneRun table (because this individual has 2 alleles total)
                            ALLELE_ONE_RUN_DATA(t,l) += 2;
                        }
                    }
                    // if this individual is female:
                    if (tree[i].sex == 2)
                    {
                        // for every cohort:
                        for (j = 0; j < tree[i].cohortNumber; j++)
                        {
                            // calculate how many cohorts after the first this is
                            t = tree[i].cohort[j] - startingCohort;
                            // add 2 to the entry for this cohort and this allele in the alleleOneRun table (because this individual has only 1 allele total)
                            ALLELE_ONE_RUN_DATA(t,l) += 1;
                        }
                    }
                }
            }
        }
    }
}

// calc_freq_one_run_data_XY function to calculate the number of alleles for each cohort for sex chromosome loci in an XY system
// This is needed to calculate allele frequency for data
void calc_freq_one_run_data_XY(int startingCohort, int * print_filter) {
    // initialize variables
    int i, j, k, l, t;

    // for every individual in the tree:
    for (i = 0; i < treeLength; i++)
    {
        // if we're including this individual
        if (print_filter[i] == 1)
        {
            // for every allele at this locus
            for (l = 0; l < nAlleles; l++)
            {
                // number of copies of this allele = 0
                k = 0;
                // if we have data from the year this individual was born:
                if (tree[i].cohortNumber > 0)
                {
                    // if this individual is female:
                    if (tree[i].sex == 1)
                    {
                        // for every cohort:
                        for (j = 0; j < tree[i].cohortNumber; j++)
                        {
                            // calculate how many cohorts after the first this is
                            t = tree[i].cohort[j] - startingCohort;
                            // add 2 to the entry for this cohort and this allele in the alleleOneRun table (because this individual has 2 alleles total)
                            ALLELE_ONE_RUN_DATA(t,l) += 2;
                        }
                    }
                    // if this individual is male:
                    if (tree[i].sex == 2)
                    {
                        // for every cohort:
                        for (j = 0; j < tree[i].cohortNumber; j++)
                        {
                            // calculate how many cohorts after the first this is
                            t = tree[i].cohort[j] - startingCohort;
                            // add 2 to the entry for this cohort and this allele in the alleleOneRun table (because this individual has only 1 allele total)
                            ALLELE_ONE_RUN_DATA(t,l) += 1;
                        }
                    }
                }
            }
        }
    }
}

void calc_freq() {
    int a, i, j;
    for (i = 0; i < nCohorts; i++) {
        if (cohortIndicator[i] == 1) {
            for (j = 0; j < nAlleles; j++) {
                a = ALLELE_ONE_RUN(i,j);
                if (a) {
                    ALLELE(i,j,a) = ALLELE(i,j,a) + 1;
                }
            }
        }
    }
}


void show_tree(struct ind *tree) {
    int i,j;
    for (i = 0; i < treeLength; i++) {
        printf("%d %d %d %d ", tree[i].i, *tree[i].a1, *tree[i].a2, tree[i].cohortNumber);
        for (j = 0; j < tree[i].cohortNumber; j++) {
            printf("%d ", tree[i].cohort[j]);
        }
        printf("\n");
    }
}

int main(int argc, char** argv) {

    FILE *fp;
    char line[1023];
    char *pt;
    int existingCohortsNumber, existingCohort, i, j, k, startingCohort, nDrops, print_long_format, allele_to_print_in_long_format, print_filter_switch, chrom_type;
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_taus2;
    gsl_rng *rnd = gsl_rng_alloc(T);
    allele_to_print_in_long_format = 2; // allele 2
    print_long_format = 0;

    /* Initialize the GSL generator with time */
    gsl_rng_set(rnd, time(NULL)); 

   /* handling first line */
    print_filter_switch = 1;
    fp = fopen(argv[1], "r");
    int cur_line = 0;
    fgets(line, 1023, fp);
    nCohorts = atoi(strtok(line,","));
    nAlleles = atoi(strtok(NULL,","));
    nIndividuals = atoi(strtok(NULL,","));
    if (strcmp(argv[2], "0") == 0) {
        print_filter_switch = 0;
    }
    if ( argc == 6 && (strcmp(argv[3], "D") == 0 || strcmp(argv[3], "d") == 0) )
    {
        nDrops = 0;
    }
    else if ( argc == 6 && (strcmp(argv[3], "R") == 0 || strcmp(argv[3], "r") == 0) )
    {
        nDrops = atoi(strtok(NULL,","));
    } 
    if (argc == 6 && (strcmp(argv[4], "L") == 0 || strcmp(argv[4], "l") == 0))
    {
        print_long_format = 1;
    } 
    else if (argc == 6 && (strcmp(argv[4], "S") == 0 || strcmp(argv[4], "s") == 0))
    {
        print_long_format = 0;
    }
    if (argc == 6 && (strcmp(argv[5], "A") == 0 || strcmp(argv[5], "a") == 0))
    {
        // set chrom_type to 0 for autosomes
        chrom_type = 0;
    } 
    else if (argc == 6 && (strcmp(argv[5], "Z") == 0 || strcmp(argv[5], "z") == 0))
    {
        // set chrom_type to 1 for ZW sex chromosomes
        chrom_type = 1;
    }
    else if (argc == 6 && (strcmp(argv[5], "X") == 0 || strcmp(argv[5], "x") == 0))
    {
        // set chrom_type to 2 for XY sex chromosomes
        chrom_type = 2;
    }
    treeLength = nIndividuals;

    cohortIndicator = malloc(nCohorts * sizeof(int));
    for (i = 0; i < nCohorts; i++) { cohortIndicator[i] = 0;}
    
    fgets(line, 1023, fp);
    existingCohortsNumber = atoi(strtok(line,","));
    for( i = 0; i < existingCohortsNumber; i++)
    {
        existingCohort = atoi(strtok(NULL,","));
        if (i == 0) { startingCohort = existingCohort; }
        j = existingCohort - startingCohort;
        cohortIndicator[j] = 1;
    }
    print_filter = malloc (nIndividuals * sizeof(int));
    a = malloc (nAlleles * sizeof(int));
    for (i = 0; i < nAlleles; i++) { a[i] = i;}
    allele = malloc ( nCohorts * nAlleles * (nIndividuals + 1) * sizeof(int));
    alleleOneRun = malloc (nCohorts * nAlleles * sizeof(int));
    alleleOneRunData = malloc (nCohorts * nAlleles * sizeof(int));
    tree = malloc (sizeof (struct ind) * treeLength );

//    allele[o3d(1,1,1)] = 1;
    
    while(cur_line < treeLength) {
        fgets(line, 1023, fp);
        fill_tree(line, print_filter, print_filter_switch);
        cur_line += 1;
   }
    fclose(fp);
    init_freq();

//  header line
    if(nDrops == 0) {
 	    printf("cohort allele allele_count all_alleles_count frequency_of_allele cohort_year\n");
    } 
    else {
    	if (print_long_format == 0) {
   	        printf("cohort allele allele_count all_alleles_count\n");
    	}
    	else {
            printf("run ");
            for (i = 0; i < nCohorts; i++) {
                printf("%d ",i);
            }
            printf("\n");    	
    	}
    }

/* when nDrops == 0, we calculate the allele counts of the original data */
    if (nDrops == 0)
    {
        init_freq_one_run();
        if (chrom_type == 0)
        {
            calc_freq_one_run_auto(startingCohort, print_filter);
            calc_freq_one_run_data_auto(startingCohort, print_filter);
        }
        else if (chrom_type == 1)
        {
            calc_freq_one_run_ZW(startingCohort, print_filter);
            calc_freq_one_run_data_ZW(startingCohort, print_filter);
        }
        else if (chrom_type == 2)
        {
            calc_freq_one_run_XY(startingCohort, print_filter);
            calc_freq_one_run_data_XY(startingCohort, print_filter);
        }
        calc_freq();
        for (i = 0; i < nCohorts; i++) {
            for (j = 0; j < nAlleles; j++) {
                for (k = 0; k < treeLength + 1; k++) {
                    if (ALLELE(i,j,k)) {
                        printf("%d %d %d %d %f %d\n", i, j, k, ALLELE_ONE_RUN_DATA(i,j),( 0.0 + k)/ALLELE_ONE_RUN_DATA(i,j), i + startingCohort );
                    }
                }
            }
        }
    }
    for (i = 0; i < nDrops; i++) {
        init_freq_one_run();
        if (chrom_type == 0)
        {
            drop_genes_faster_auto(rnd);
            calc_freq_one_run_auto(startingCohort, print_filter);
        }
        else if (chrom_type == 1)
        {
            drop_genes_faster_ZW(rnd);
            calc_freq_one_run_ZW(startingCohort, print_filter);
        }
        else if (chrom_type == 2)
        {
            drop_genes_faster_XY(rnd);
            calc_freq_one_run_XY(startingCohort, print_filter);
        }
        calc_freq();
        if (print_long_format == 1)
        {
            print_freq_one_run_v2(i, allele_to_print_in_long_format);
        }
    }
    if (print_long_format == 0 && nDrops > 0)
    {
        for (i = 0; i < nCohorts; i++) {
            for (j = 0; j < nAlleles; j++) {
                for (k = 0; k < treeLength + 1; k++) {
                    if (ALLELE(i,j,k)) {
                        printf("%d %d %d %d\n", i, j, k, ALLELE(i,j,k));
                    }
                }
            }
        }
    }
    return(1);
}
