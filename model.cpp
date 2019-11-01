//
// Created by verniy on 18-10-29.
//

#include "model.h"

namespace ntoolkit{
  Neuron_net_base::Neuron_net_base()
          : neuron({}), net(nullptr){
  }

  Neuron_net_base::~Neuron_net_base(){
      Clear();
  }

  std::vector<double> Neuron_net_base::Refresh()
  {
      std::vector<double> result;
      if (neuron.empty())
          return result;
      unsigned long neuron_num = neuron.size();
      auto * other_neuron = new double[neuron_num]{};
      if (net != nullptr) {
          tbb::parallel_for(0, (int) neuron_num, 1, [&](int i) {
              for (int j = 0; j < neuron_num; j++) {
                  other_neuron[i] += net->Weight(std::make_pair(i, j)) * neuron.at(j)->Measure();
              }
          });
      }

      // should sum other neurons' spkie-weight before refresh each neuron
      result.resize(neuron_num);
      tbb::parallel_for(0, (int)neuron_num, 1, [&](int i) {
          neuron[i]->Refresh({other_neuron[i],0.0});
          result[i] = neuron.at((unsigned long)i)->Measure();
      });

      for (int i = 0; i < neuron_num; i++){
          neuron[i]->Refresh({other_neuron[i]});
          result.push_back(neuron.at(i)->Measure());
      }
      delete[] other_neuron;
      return result;
  }

  std::vector<double> Neuron_net_base::Refresh(const std::vector<double> &ext) {
      std::vector<double> result;
      if (neuron.empty())
          return result;
      unsigned long neuron_num = neuron.size();
      auto * other_neuron = new double[neuron_num]();
      if (net != nullptr) {
          tbb::parallel_for(0, (int) neuron_num, [&](int i) {
              for (int j = 0; j < neuron_num; j++) {
                  other_neuron[i] += net->Weight(std::make_pair(i, j)) * neuron.at(j)->Measure();
              }
          });
      }
      // should sum other neurons' spkie-weight before refresh each neuron
      result.resize(neuron_num);
      tbb::parallel_for(0, (int)neuron_num, 1, [&](int i) {
          neuron[i]->Refresh({other_neuron[i], ext.at((unsigned long)i)});
          result[i] = neuron.at((unsigned long)i)->Measure();
      });
      for (int i = 0; i < neuron_num; i++){
          neuron[i]->Refresh({other_neuron[i]});
          result.push_back(neuron.at(i)->Measure());
      }
      delete[] other_neuron;
      return result;
  }

  bool Neuron_net_base::GenNetReport(std::string filename) {
      return true;
  }

  void Neuron_net_base::Clear(){
      for( auto &it : neuron){
          delete it;
          it = nullptr;
      }
      neuron.clear();
      delete net;
      net = nullptr;
  }

/*
AiharaBrainNet::AiharaBrainNet()
        :Neuron_net_base()
{
}

void AiharaBrainNet::InitialNeuron(const vector<double>& paras_name_list, const vector<double>& initial)
{
    auto *para = new double[4]{ paras_name_list.at(0), paras_name_list.at(1), paras_name_list.at(2), paras_name_list.at(3) };
    //para[3] = paras_name_list.at(3);
    //para[2] = paras_name_list.at(2);
    //para[1] = paras_name_list.at(1);
    //para[0] = paras_name_list.at(0),;
    for (const double& i : initial)
    {
        Neuron* temp;
        temp = new Aihara(i,para,paras_name_list.size());
        neuron.push_back(temp);
    }
    delete[]para;
}

void AiharaBrainNet::InitialNet(const double * weight)
{
    net = new Brain_net(neuron.size(), weight);
}

HHBrainNet::HHBrainNet()
        :Neuron_net_base(6)
{
}


void HHBrainNet::InitialNeuron(const vector<double>& paras_name_list, const vector<double>& initial)
{
    for (const double& i : initial)
    {
        Neuron* temp;
        temp = new Hodgkin_huxley(i, paras_name_list.at(0), paras_name_list.at(1), paras_name_list.at(2), paras_name_list.at(3), paras_name_list.at(4), paras_name_list.at(5));
        neuron.push_back(temp);
    }
}

void HHBrainNet::InitialNet(const double * weight)
{
    net = new Brain_net(neuron.size(), weight);
}


FHNBrainNet::FHNBrainNet()
        :Neuron_net_base(5)
{
}

void FHNBrainNet::InitialNeuron(const vector<double>& paras_name_list, const vector<double>& initial)
{
    for (int i=0; i<initial.size(); i++)
    {
        Neuron* temp;
        temp = new Fitzhugh_nagumo(paras_name_list.at(0), paras_name_list.at(1), paras_name_list.at(2), paras_name_list.at(3), paras_name_list.at(4));
        neuron.push_back(temp);
    }
}

void FHNBrainNet::InitialNet(const double * weight)
{
    net = new Brain_net(neuron.size(), weight);
}

FHNBrainNet2::FHNBrainNet2()
        :Neuron_net_base(3)
{
}

void FHNBrainNet2::InitialNeuron(const vector<double>& paras_name_list, const vector<double>& initial)
{
    for (int i=0; i<initial.size(); i++)
    {
        Neuron* temp;
        temp = new Fitzhugh_nagumo2(0.6, 0.0, 0.1);
        neuron.push_back(temp);
    }
}

void FHNBrainNet2::InitialNet(const double * weight)
{
    net = new Brain_net(neuron.size(), weight);
}
*/
}
