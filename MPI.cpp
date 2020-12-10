/*
Ciprian Alex Cornila
ECE 4122 Advanced Programming Techniques
11/09/2019
*/
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <stdio.h>
#include "mpi.h"
#include <stddef.h>
#include <limits>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <random>
#include <omp.h>
#define NUM_THREADS 4
#define  MASTER		0
using namespace std;

//struct to store every YJ and buzzy data
typedef struct ship
{
  int tag;
  int status;
  double pos[3];
  double vel[3];
  double accel[3];
} ship;

//global function helpers
double randF();
bool NextToDock(int rank, ship a[]);
void PrintData(ship ship[]);
double Vmag(ship ship[], int r);
bool dockedYJ(bool YJstatus);




int main(int argc, char**argv)
{
  int  numtasks, rank, rc, time;
  ship shipArray[8];
  double maxForce;
  bool dockedYJ = false;
  bool nextYJ = false;
  bool dockPos = false;
  string dat;
  vector<string> Vec;

  //read in data filename replace this
  std::string filename = "in.dat";
  //initilizing MPI
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
  {
    printf("MPI not initialized\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //to send data via MPI we need to serialize structs by definging its layout to MPI
  const int nitems=5; // we need 5 items, members of struct
  int blocklengths[5] = {1, 1, 3, 3, 3}; //size of each item to send
  MPI_Datatype types[5] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}; //type of the data MPI to handle
  MPI_Datatype mpi_ship_type;//type of struct to repack the data
  MPI_Aint offsets[5];
  offsets[0] = offsetof(ship, tag);
  offsets[1] = offsetof(ship, status);
  offsets[2] = offsetof(ship, pos);
  offsets[3] = offsetof(ship, vel);
  offsets[4] = offsetof(ship, accel);
  MPI_Type_struct(nitems, blocklengths, offsets, types, &mpi_ship_type);
  MPI_Type_commit(&mpi_ship_type);
 //https://stackoverflow.com/questions/9864510/struct-serialization-in-c-and-transfer-over-mpi


   //reading data from input file
   std::string line;
   ifstream file (filename);
   file >> dat;
   time = stoi(dat);
   file >> dat;
   maxForce = stof(dat);

  //Master process initializes the array of shipArray
  if (rank == MASTER)
  {
    //read in data for master process
    std::string line;
    ifstream file (filename);
    int counter = 1;
    if (file.is_open())
    {
      while ( getline (file,line) )
      {
        //read the game time
        if (counter == 1)
        {
          time = stoi(line);
        }
        //read the max force
        else if (counter == 2)
        {

          maxForce = stof(line);
        }
        // lines 3+: direction and speed of each ship
        else
        {
          //store data in Vec
          Vec.clear();
          string word = "";

          //parallel region per process starts here
          //#pragma omp parallel for
          for (auto x : line)
          {
            if (x == ' ')
            {
              Vec.push_back(word);
              word = "";
            }
            else
            {
              word = word + x;
            }
          }
          //https://www.geeksforgeeks.org/split-a-sentence-into-words-in-cpp/
          //#pragma omp critical  //end of parallel region
          //skip spaces and push all words
          Vec.push_back(word);
          int YJcounter = counter -3;

          //read in data and store into shipArray
          for (int i=0; i<3; i++)
          {
            //set all active initially
            shipArray[YJcounter].status = 1;
            shipArray[YJcounter].accel[i] = 0;
            //initial pos from file
            shipArray[YJcounter].pos[i] = stof(Vec[i]);
            //initial vel from file
            shipArray[YJcounter].vel[i] = stof(Vec[3]) * stof(Vec[4+i]);
            shipArray[YJcounter].tag = YJcounter;
          }
        }
        counter++;
      }
      file.close();
    }
    else cout << "err file";
  }


 //parallel region for processes starts here
  #pragma omp parallel
  //looping over the  entire time allocated
  while(--time)
  {
    //master thread
    if (rank == MASTER)
    {
      MPI_Status status;

      //loop over shipArray & send/recv data via MPI communicator
      for (int i=1; i < 8; i++)
      {
        MPI_Send(&shipArray, 8, mpi_ship_type, i, 0, MPI_COMM_WORLD);

        MPI_Recv(&shipArray[i], 1, mpi_ship_type, i, 0, MPI_COMM_WORLD, &status);
      }

      //move BattlestarBuzzy by adjusting x,y,z
      for (int i=0; i<3; i++)
      {
        //pos = initial_pos + vel * t
        shipArray[0].pos[i] = shipArray[0].pos[i] + shipArray[0].vel[i];
      }

       //print shipArray current stats on every itteration
       PrintData(shipArray);
       //std::cout << time<<":------------------------------------"<<std::endl;

      //check for non docked YJ
      int dockedYJ = {0};
      for (int i=1; i<8; i++)
      {
        if (shipArray[i].status != 1)
          dockedYJ ++;
      }
      //break out of the loop once all YJ are docked
      if (dockedYJ == 7)
      {
        exit(MPI_SUCCESS);//exit MPI if all YJ are docked
      }
    }


    //YJ parameters
    if(rank != MASTER)
    {
      MPI_Status status;
      //rec from main process
      MPI_Recv(&shipArray, 8, mpi_ship_type, 0, 0, MPI_COMM_WORLD, &status);

      //store distances from the YJ and calc magnitude
      double dist[8];
      for (int i=0; i<8; i++)
      {
        //distance = magnitude of the squared diffreces
        dist[i] = {0};
        for (int j=0; j<3; j++)
        { //loop over the coordinates x, y, z and get dist Xyj - Xbuzzy
          dist[i] += pow((shipArray[i].pos[j] - shipArray[rank].pos[j]),2);
        }
        dist[i] = sqrt(dist[i]);//mag of dist
      }

      //check for collision conditions between YJ
      for (int i=1; i<8; i++)
      { //when distance is < 250 and both active -> destroyed
        if (dist[i] < 250 && i != rank && shipArray[i].status == shipArray[rank].status)
        {
          shipArray[rank].status = 0;
        }
      }

      //checking if docking conditions are met
      if (dist[0] < 50)
      {
        double theta;
        double VelYJ;
        double VelBuzz;
        double dotProd;

        for (int i=0; i<3; i++)
        {
          VelBuzz += pow(shipArray[0].vel[i],2);
          VelYJ += pow(shipArray[rank].vel[i],2);
        }
        //mag = sqrt(mag)
        VelBuzz = sqrt(VelBuzz);
        VelYJ = sqrt(VelYJ);

        dotProd = {0};
        for (int i=0; i<3; i++)
        {
          dotProd += shipArray[0].vel[i] * shipArray[rank].vel[i];
        }
        //det angle of orientation for docking
        theta = dotProd / (VelBuzz * VelYJ);

        // if theta < 0.8 & VelYJ < 1.1*VelBuzz ready to dock
        //else YJ crushed
        if (theta > 0.8 && VelYJ < 1.1 * VelBuzz)
        {
          shipArray[rank].status = 2;
        }
        else
        {
          shipArray[rank].status = 0;
        }
      }

      //get next YJ to dock and allign it with buzzy
      if (shipArray[rank].status == 1)
      {
        //get nextYJ to dock use nextYJ flag to  only work with current YJ
        //that is active and ignore other active YJ
        nextYJ = true;
        for (int i=1; i<rank; i++)
        {
          if (shipArray[i].status == 1)
          {
            nextYJ = false;
          }
        }
        //allign YJ in pos to dock
        if (nextYJ && dockPos == false)
        {
          shipArray[0].pos[2] -= 500;
          for (int i=0; i<3; i++)
          {
            //adjust for to match Buzzy
            double desiredForce = (shipArray[0].pos[i]-shipArray[rank].pos[i]-shipArray[rank].vel[i]);
            shipArray[rank].accel[i] = desiredForce*1.6;
          }
          //check if YJ pos matches buzzy
          //determine if the xyz pos components of the docking distance are met
          double xDistance = abs(shipArray[0].pos[0] - shipArray[rank].pos[0]);
          double yDistance = abs(shipArray[0].pos[1] - shipArray[rank].pos[1]);
          double zDistance = abs(shipArray[0].pos[2] - shipArray[rank].pos[2]);
            //check relative pos to buzzy based on precalculated values
            if (xDistance < 3.7 && yDistance < 3.7 && zDistance < 250)
          {
            dockPos = true;
          }
        }
        //next YJ ready to dock
        else if(nextYJ && dockPos == true)
        {
          //det magnitude vel sqrt(sum of squared vec_components) 3Dpythagorean
          double VelBuzz = {0};
          double VelYJ = {0};
          for (int i=0; i<3; i++)
          {
            VelBuzz += pow(shipArray[0].vel[i],2);
            VelYJ += pow(shipArray[rank].vel[i],2);
          }
          VelBuzz = sqrt(VelBuzz);
          VelYJ = sqrt(VelYJ);

          //calc and apply desired force to match buzzy
          for (int i=0; i<3; i++)
          {
            double desiredForce;
            //motion eq to calc desired force
            if (i == 2)
            {
              desiredForce = shipArray[0].vel[i]*1.08 - shipArray[rank].vel[i];
            }
            else
            {
              desiredForce = shipArray[0].vel[i] - shipArray[rank].vel[i];
            } //check for overpowering force and choose the smallest of the two
            if(desiredForce < maxForce)
              {
                shipArray[rank].accel[i] = desiredForce;
              }
              else
              {
                shipArray[rank].accel[i] = maxForce;
              }
          }
        }

        else
        {
          //keep not currently docking YJ behind
          shipArray[0].pos[2] -= 900;// anything less that 700 crushes some YJ
          //adjust docking YJ Force to match Buzzy
          double desiredForce = (shipArray[0].pos[2]-shipArray[rank].pos[2]-shipArray[rank].vel[2]);
          if(desiredForce < maxForce)
          {
            shipArray[rank].accel[2] = desiredForce;
          }
          else
          {
            shipArray[rank].accel[2] = maxForce;
          }
            shipArray[rank].accel[1] = 0;
            shipArray[rank].accel[0] = 0;
        }
        //apply random force to current accel components
        for (int i=0; i<3; i++)
        {
          shipArray[rank].accel[i] *= randF();
          //new pos and accel based on rand force gen
          shipArray[rank].pos[i] += shipArray[rank].vel[i];
          shipArray[rank].vel[i] += shipArray[rank].accel[i];
          shipArray[rank].pos[i] += .5 * shipArray[rank].accel[i];
        }
      }
       //send to master
       MPI_Send(&shipArray[rank], 1, mpi_ship_type, 0, 0, MPI_COMM_WORLD);
    }
  }
  #pragma omp critical
  MPI_Finalize();
  exit(1);
}//end main



//random no generator that simulates damage thrust
double randF()
{
  std::random_device rd; //hardware gen rand
  std::mt19937 eng(rd());
  std::uniform_int_distribution<> distr(80, 120);
  double ran = static_cast<double>(distr(eng))/100;
  return ran;
}

//get next YJ to dock
bool NextToDock(int rank, ship a[])
{
  bool next = true;
  for (int i=1; i<rank; i++)
  {
    if (a[i].status == 1)
    {
      bool next = false;

    }
  }
  return next;
}

// computes magnitude of xyz  vec components
double Vmag(ship ship[], int r)
{
  double Mag;
  for (int i=0; i<3; i++)
  {
    Mag += pow(ship[r].vel[i],2);

  }
  //mag = sqrt(mag)
   return Mag = sqrt(Mag);

}
//prints shipArray data
void PrintData(ship ship[])
{

  for (int i=1; i<8; i++)
  {
    cout << setprecision(6) << ship[i].tag << "," << ship[i].status;
    for (int j=0; j<3; j++)
    {
      cout << "," <<setprecision(6) << ship[i].pos[j];
      //cout << "shipArray" << std::endl;
    }
    for (int j=0; j<3; j++)
    {
      cout << "," << ship[i].accel[j]*10000;
    }
    cout << endl;
  }
}
