//
// Created by verniy on 18-10-29.
//
/*!
 * @file
 * @author Tooko Madobe <madobet@outlook.com>
 * @version 0.1
 *
 * @section LICENSE
 *
 * @section DESCRIPTION
 *
 * 单一神经元模型
 */


#include "neuron.h"

namespace ntoolkit{
  /*!
   *
   * @param subclass_para_list
   * @param para_inputs
   */
  Neuron::Neuron(
          // 神经元基类
          const std::vector<const std::string> &subclass_para_list,
          // 子类所需的参数项组成的向量，每个元素是一个参数名
          std::initializer_list<double> para_inputs
          // 各参数名对应的参数值，可以在最后带上初值，如果不指定，则取0值初始
  )
          :paras_name_list(subclass_para_list)
  {
      if(para_inputs.size() < paras_name_list.size())     throw "模型参数太少";
      if(para_inputs.size() > (paras_name_list.size()+1)) throw "模型参数太多";
      auto *p = para_inputs.begin();
      for(int i=0; i<paras_name_list.size(); i++){
          parameters[paras_name_list.at(i)] = *p; p++;
      }
      // 如果不指定，则取0值初始
      if(p==para_inputs.end())    spike = 0;
      else                        spike = *p;
  }

  std::string Neuron::PrintUsage() const {
      std::string res("");
      for(const std::string& it : paras_name_list) res += (it + ' ');
      res += "initial-status";
      return res;
  }
  std::string Neuron::PrintParaList() const {
      std::string res("");
      for(const std::string& it : paras_name_list) res += (it + '_');
      return res;
  }
  Aihara::Aihara(std::initializer_list<double> paras_un_init)
          : Neuron({"kf", "kr", "ai", "alpha" }, paras_un_init), eta_i(0), zeta_i(0), stat_inter(0){
  }
  Aihara::Aihara(const Aihara& src): eta_i(src.eta_i), zeta_i(src.zeta_i), stat_inter(src.stat_inter){
      spike = src.spike;
      parameters = src.parameters;
  }
  Aihara& Aihara::operator=(const Aihara& src){
      this->parameters = src.parameters;
      this->eta_i = src.eta_i;
      this->zeta_i = src.zeta_i;
      this->stat_inter = src.stat_inter;
      return *this;
  }
  Aihara& Aihara::operator*(double dev){
      eta_i += eta_i*dev;
      zeta_i += zeta_i*dev;
      stat_inter = eta_i + zeta_i;
      spike = math::Sigmoid(stat_inter);
      return *this;
  }
  void Aihara::Refresh(std::initializer_list<double> input){
      static const double kf(parameters["kf"]), kr(parameters["kr"]), alpha(parameters["alpha"]), ai(parameters["ai"]);
      eta_i = kf*eta_i + std::accumulate(input.begin(), input.end(), (double)0.0);
      zeta_i = kr*zeta_i - alpha*spike + ai;
      stat_inter = eta_i + zeta_i;
      spike = math::Sigmoid(stat_inter);
  }

  const double Hodgkin_huxley::C_m = 1;
  const double Hodgkin_huxley::g_Na_coff = 120;
  const double Hodgkin_huxley::g_K_coff = 36;
  const double Hodgkin_huxley::g_L = 0.3;
  const double Hodgkin_huxley::E_rest = -60;
  const double Hodgkin_huxley::rho_Na = 60;
  const double Hodgkin_huxley::rho_K = 18;
  Hodgkin_huxley::Hodgkin_huxley(std::initializer_list<double> paras_und_init)
          : Neuron({"ion_Na_i", "ion_Na_e", "ion_K_i", "ion_K_e", "E_L", "dt"}, paras_und_init){
      dV = spike;
      spike += E_rest;
      // there is a bug  with Neuron class initialize
      E_Na = 26 * log(parameters["ion_Na_e"]/parameters["ion_Na_i"]);
      E_K  = 26 * log(parameters["ion_K_e"] /parameters["ion_K_i"] );
      V_Na = E_rest - E_Na;
      V_K  = E_rest - E_K;
      V_L  = E_rest - parameters["E_L"];
      m = alpha_m(dV) / (alpha_m(dV) + beta_m(dV));
      // V_0 = dV = 0 mV than m ~= 0.052932
      h = alpha_h(dV) / (alpha_h(dV) + beta_h(dV));
      // V_0 = dV = 0 mV than h ~= 0.59612
      n = alpha_n(dV) / (alpha_n(dV) + beta_n(dV));
      // V_0 = dV = 0 mV than n ~= 0.31768
  }
  void Hodgkin_huxley::Refresh(std::initializer_list<double> input){
      static const double E_L(parameters["E_L"]), dt(parameters["dt"]);
      double x = spike;
      double f = (- g_Na_coff * pow(m, 3)*h*(x - E_Na) \
                - g_L * (x - E_L) \
                - g_K_coff * pow(n, 4) * (x - E_K) + std::accumulate(input.begin(), input.end(), (double)0.0)) / C_m;
      double k_1 = f;
      double k_2 = f + dt / 2 * k_1;
      double k_3 = f + dt / 2 * k_2;
      double k_4 = f + dt * k_3;
      spike = x + dt / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
      f = alpha_m(-dV) * (1 - m) - beta_m(-dV)*m;
      k_1 = f;
      k_2 = f + dt / 2 * k_1;
      k_3 = f + dt / 2 * k_2;
      k_4 = f + dt * k_3;
      m = m + dt / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
      f = alpha_h(-dV) * (1 - h) - beta_h(-dV)*h;
      k_1 = f;
      k_2 = f + dt / 2 * k_1;
      k_3 = f + dt / 2 * k_2;
      k_4 = f + dt * k_3;
      h = h + dt / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
      f = alpha_n(-dV) * (1 - n) - beta_n(-dV)*n;
      k_1 = f;
      k_2 = f + dt / 2 * k_1;
      k_3 = f + dt / 2 * k_2;
      k_4 = f + dt * k_3;
      n = n + dt / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
      dV = spike - E_rest;
      /*
      V_Na = E_Na - V_rest;
      V_K = E_K - V_rest;
      V_L = E_L - V_rest;
      */
  }
  
  Fitzhugh_nagumo::Fitzhugh_nagumo(std::initializer_list<double> paras_und_init)
          : Neuron({"I", "tau", "a", "b", "dt"}, paras_und_init){
      static const double a(parameters["a"]), b(parameters["b"]), I(parameters["I"]);
      const double err = 1e-14;
      double va = b;      //                      b         v^3
      double vb = 0;      //                      0         v^2
      double vc=3*(1-b);  //                      3(1-b)    v
      double vd=3*(a-b*I);//          +           3(a-bI)   1
      //--------------------------------------------------------
      // b   v^3 + 0      v^2 + 3(1-b)        v + 3(a-bI)    = 0
      double wa = pow(b,3);
      double wb = a*pow(b,2);
      double wc = pow(a,2)*b+3*(1-b);
      double wd = 3*(a-I)-pow(a,3);
      //--------------------------------------------------------
      // b^3 w^3 + a(b^2) w^2 + (a^2)b+3(1-b) w + 3(a-I)-a^3 = 0
      spike = math::NewtonCubicEquRT(-0.5,err,va,vb,vc,vd);
      w = math::NewtonCubicEquRT(-0.5,err,wa,wb,wc,wd);
  }
  void Fitzhugh_nagumo::Refresh(std::initializer_list<double> input)
  {
      static const double a(parameters["a"]), b(parameters["b"]), tau(parameters["tau"]), dt(parameters["dt"]);
      double x = spike;
      double f = x - pow(x,3) - w + std::accumulate(input.begin(), input.end(), (double)0.0);
      double k_1 = f;
      double k_2 = f + dt / 2 * k_1;
      double k_3 = f + dt / 2 * k_2;
      double k_4 = f + dt * k_3;
      spike = x + dt / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
      //
      f = (x + a - b*w) / tau;
      k_1 = f;
      k_2 = f + dt / 2 * k_1;
      k_3 = f + dt / 2 * k_2;
      k_4 = f + dt * k_3;
      w = w + dt / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
  }
  Fitzhugh_nagumo2::Fitzhugh_nagumo2(std::initializer_list<double> paras_und_init)
          : Neuron({"delta", "dt"}, paras_und_init), epsilon(1e-6){
      w = spike-pow(spike,3)/3;
  }
  void Fitzhugh_nagumo2::Refresh(std::initializer_list<double> input){
      static const double dt(parameters["dt"]);
      double x = spike;
      double f = (x - pow(x,3) - w + std::accumulate(input.begin(), input.end(), (double)0.0))/epsilon;
      double k_1 = f;
      double k_2 = f + dt / 2 * k_1;
      double k_3 = f + dt / 2 * k_2;
      double k_4 = f + dt * k_3;
      spike = x + dt / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
      f = spike;
      k_1 = f;
      k_2 = f + dt / 2 * k_1;
      k_3 = f + dt / 2 * k_2;
      k_4 = f + dt * k_3;
      w = w + dt / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
  }

}   // namespace ntoolkit
