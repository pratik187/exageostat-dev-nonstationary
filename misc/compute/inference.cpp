#include <iostream>
#include "cppflow/cppflow.h"
#include "../include/inference.hpp"
// #include "../include/MLE_misc.h"
#include <bits/stdc++.h>
#include <unistd.h>
#include <stdio.h>
#include <limits.h>
using namespace std;

pair<int, int> findGrid(double x, double y, double startx, double endx, double starty, double endy)
{
    // n=3600 ONLY. 
    // cout << "I am here inside findGrid \n "; 
    // printf("%f__%f__%f__%f__%f__%f\n",x,y,startx,endx,starty,endy);
    int gridi = floor(((y - starty) * 100) / (endy - starty));
    int gridj = floor(((x - startx) * 100) / (endx - startx));

    pair<int, int> ans = {min(99, gridi), min(99, gridj)};
    return ans;
}

//
// startx, starty, endx, endy denotes the square subregion that needs to be evaluated
// Function returns a normalised image of size 60 x 60 with each pixel containing average of all locations lying inside that pixel.

vector<double> preprocess(double *x, double *y, double* Z_values, double startx, double endx, double starty, double endy, double startz, double endz , bool z_user_range , int N)
{
    // cout << "I am here \n "; 
    int sz = N;
    // cout << "I am here \n "; 
    vector<pair<double, int>> grids(10000, {0.0, 0});
    // cout << "I am here \n "; 
    //  X--
    //  ---
    //  ---
    //  Number of locations in the (i,j)th subgrid = grids[i*60 + j].second
    //  Sum of Z-values of all locations in (i,j)th subgrid = grids[i*60 + j].first
    printf("%d\n",sz);
    for (int i = 0; i < sz; i++)
    {
        pair<int, int> subgrid = findGrid(x[i], y[i], startx, endx, starty, endy);
        int subgridi = subgrid.first;
        int subgridj = subgrid.second;
        grids[subgridi * 100 + subgridj].second++;
        grids[subgridi * 100 + subgridj].first += Z_values[i];
    }
    
    // char * nFileLoc1  = (char *) malloc(100 * sizeof(char));
    // snprintf(nFileLoc1, 100, "/home/nagp/processed_data1");
    // FILE *fptr2;
    // fptr2 = fopen(nFileLoc1,"w");
    // for(int i=0; i<10000; i++){
    //     fprintf(fptr2,"%f,%d\n",grids[i].first,grids[i].second);
    // }
    // fclose(fptr2);
    
    // double mean_Z = accumulate(Z_values[0], Z_values[N-1], 0.0) / sz; 
    // double sum_z = 0.0;
    // for (int i = 0;i < N; i++){
    //     sum_z = sum_z + Z_values[i];
    // }
    // double mean_Z = sum_z/sz;
    // cout << "Mean value of Z:" << mean_Z << endl; 
    
    vector<double> processed_grid(10000);
    
    for(int i = 0; i < 100; i++){
        for (int j= 0; j < 100; j++)
        {
            int count=0;
            double sum=0;
            if (grids[i * 100 + j].second == 0){
                // cout << "I am here \n "; 
                // cout << "value of i : " << i << "value of j : " << j << endl;
                for(int k=i-1; k<=i+1;k++){
                    for(int l=j-1;l<=j+1;l++){
                        // cout << "value of k : " << k << "value of l : " << l << endl;
                        if(k>=0 && k<100 && l>=0 && l<100){
                            // cout << "I am here \n "; 
                            if (grids[k * 100 + l].second != 0){
                                count++;
                                sum += grids[k * 100 + l].first/grids[k * 100 + l].second ;
                                
                            }
                            // cout << "value of k : " << k << "value of l : " << l << endl;
                        } 
                    }
                
                }
                if(count != 0){
                    processed_grid[i * 100 + j] = sum/count;
                }else{
                    for(int k=i-10; k<=i+10;k++){
                        for(int l=j-10;l<=j+10;l++){
                            // cout << "value of k : " << k << "value of l : " << l << endl;
                            if(k>=0 && k<100 && l>=0 && l<100){
                                // cout << "I am here \n "; 
                                if (grids[k * 100 + l].second != 0){
                                    count++;
                                    sum += grids[k * 100 + l].first/grids[k * 100 + l].second ;
                                }    
                            }    
                                // cout << "value of k : " << k << "value of l : " << l << endl;
                        } 
                    }
                
                }
            
            }else{
                processed_grid[i * 100 + j] = grids[i * 100 + j].first / grids[i * 100 + j].second;
            }
        }
    }
    
    // for (int i = 0; i < 3600; i++)
    // {
    //     if (grids[i].second == 0)
    //         processed_grid[i] = mean_Z;
    //     else
    //         processed_grid[i] = grids[i].first / grids[i].second;
    // }

    double grid_min = *min_element(processed_grid.begin(), processed_grid.end());
    double grid_max = *max_element(processed_grid.begin(), processed_grid.end());

    // cout << "Min element: " << grid_min <<"; Max element: " << grid_max << endl; 

    if(z_user_range == false) {
        transform(processed_grid.begin(), processed_grid.end(), processed_grid.begin(), [&grid_min, &grid_max](double c) -> double {return (c - grid_min)/(grid_max - grid_min); } ); 
    } else {
        transform(processed_grid.begin(), processed_grid.end(), processed_grid.begin(), [&startz, &endz](double c) -> double {return (c - startz)/(endz - startz); } ); 
    }
    //check transformed unit
    grid_min = *min_element(processed_grid.begin(), processed_grid.end()); 
    grid_max = *max_element(processed_grid.begin(), processed_grid.end()); 
    // cout << "Min element: " << grid_min <<"; Max element: " << grid_max << endl;

    // char * nFileLoc  = (char *) malloc(100 * sizeof(char));
    // snprintf(nFileLoc, 100, "/home/nagp/processed_data");
    // FILE *fptr1;
    // fptr1 = fopen(nFileLoc,"w");
    // for(int i=0; i<10000; i++){
    //     fprintf(fptr1,"%f\n",processed_grid[i]);
    // }
    // fclose(fptr1);
    // grids.clear();
    // grids.shrink_to_fit();
    return processed_grid;
}

int predict(double *x, double *y, double *Z_values, double startx, double endx, double starty, double endy, double startz, double endz, bool z_user_range, int N){

    
    vector<double> processed_grid = preprocess(x , y , Z_values , startx, endx, starty, endy, startz,endz, z_user_range, N);
    vector<float> processed_grid_float(processed_grid.begin(), processed_grid.end());
    
    char cwd[PATH_MAX];
    char FileLoc[100];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
    printf("Current working dir: %s\n", cwd);
    } else {
    perror("getcwd() error");
    }
    // cout << "I am here \n "; 
    snprintf(FileLoc, 100, "%s/Model_Example", cwd);
    printf("Current working dir: %s\n", FileLoc);
    cppflow::model model(FileLoc);
    // cppflow::model model("/home/nagp/exageostat-dev/misc/compute/model/SavedModel.pb");
    
    vector<int64_t> sh = {100, 100};
    auto a = cppflow::tensor( processed_grid_float , sh);
    
    auto input = cppflow::expand_dims(a, -1);
    input = cppflow::expand_dims(input, 0);
    auto output = model(input);
    auto logits = output;
    auto probabilities = cppflow::softmax(logits);
    auto prob = probabilities.get_data<float>();

    //cout<<"Probabilities : "<<probabilities<<endl;
    cout<<prob.size()<<endl;
    // cout<<"Probability of being stationary:" << prob[0] << "; being nonstationary: " << prob[1] <<endl;

    double proba = prob[1]*10000;
    
    // cout << "probability before :" << probability << endl;
    // probability = proba;
    // proba=56.76;
    //cout << typeid(proba).name() << endl;
    cout << "(inside c++)probability after :" << proba << endl;
    //int awwer1 = 5;
    
    return (int)proba;

}

// int main()
// {

//     auto arr = {1.0, 2.0, 3.0};
//     vector<double> v(3600, 0.0);
//     v[2] = 1.0;
//     vector<int64_t> sh = {60, 60};
//     auto a = cppflow::tensor(v, sh);

//     // auto b = cppflow::fill({ 3 , 3 }, {1.0});
//     cppflow::model model("./model");

//     // Get all operations in the model
//     vector<string> ops = model.get_operations();
//     for (auto xx : ops)
//     {
//         std::cout << xx << std::endl;
//     }

//     // Add batch dimension
//     auto input = cppflow::expand_dims(a, -1);
//     input = cppflow::expand_dims(input, 0);

//     std::cout << input.shape() << std::endl;
//     auto output = model(input);
//     auto logits = output;

//     // Apply softmax to logits to get probabilities
//     auto probabilities = cppflow::softmax(logits);
//     auto data = output.get_data<double>();

//     std::cout << "Model output is : " << cppflow::arg_max(output, 1) << std::endl;
//     std::cout << "Model output shape is : " << output.shape() << std::endl;
//     std::cout << "Model output raw is : " << output << std::endl;
//     std::cout << "Model output raw data is : " << data[1] << std::endl;
//     std::cout << "Model output probabilities are : " << probabilities << std::endl;

//     std::cout << a << std::endl;
//     std::cout << a.shape() << std::endl;

//     return 0;
// }
