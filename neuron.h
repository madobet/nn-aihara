//
// Created by verniy on 18-10-29.
//

#ifndef NTOOLKIT_NEURON_H
#define NTOOLKIT_NEURON_H

#include <vector>
#include <map>
#include "lib.h"

namespace ntoolkit{
/* * * * * * * * * * * *
 * Virtual base class of all neuron models
 * 继承Neuron_base类以方便的创建特定的神经元模型
 * Neuron_base是一个虚基类，因此您不能直接使用它
 * 请在子类中对Refresh方法进行重载
 * 【使用方法】
 * 1.创建您需要的神经元模型类，继承虚基类Neuron_base
 * 2.在构造函数中将该模型的可变参数名（以字符串数组的形式）、
 *  以及您的构造函数接受到的可变参数列表（一个std::initializer_list<double>类型）传递给基类
 * 3.当您需要使用某一个参数的值时，像使用一个数组一样，使用parameter["name"]即可，name就是您传递给父类的某一个可变参数的名字
 * 4.override父类纯虚函数void Refresh(std::initializer_list<double> )方法
 * 5.加入和您模型相关的其他成员变量或者成员函数
 * 6.Just using it!
 *
 * 【关于设计】
 * 为什么使用std::initializer_list<double>而不是std::initializer_list<double>&，也即引用类型？
 * 因为标准库设计std::initializer_list<double>本身就是一种引用
 *
 * vector<string> paras_name_list   存储参数名的字符串
 * map<string,double> parameters    存储参数的具体数值
 * 更简单的做法似乎是使用枚举+纯粹的vector<double>，显然这样的做法应该？能提升性能，但会增加代码的维护成本和出错概率。
 * 枚举类型不同于一般变量类型，准确的说枚举不是变量的类型，而是类型的类型，枚举表明了一种类型，因此：
 * 1.首先子类无法继承；其次就算子类继承了，但每个模型的参数名称和个数都不同，而您无法“重载”一个枚举
 * 2.如果用不限定作用域的枚举类型？
 *  1) 子类多了枚举之间免不了容易产生冲突
 *  2) 您就无法在基类中完成通用的参数载入，需要在每个子类的构造函数里实现各自的参数载入（因为避开冲突所以枚举中元素的名称不同）
 * 3.如果使用限定作用于的枚举？确实，这样您就可以在父类和子类、子类和子类之间使用相同的枚举名，
 *  但是当您在基类中使用 作用域::枚举元素名 的方式实现通用的参数载入方法时，父类只会使用自身作用域内的枚举，并不会使用子类的同名枚举
 * 4.如果将通用的操作封装成一个函数？和3一样，即使子类继承了父类的该函数，函数也只会使用父类的枚举，而非子类的同名枚举
 * 综上所述，对于本程序，正确性和泛型（泛型亦可减少出错的概率）高于运行速度，因此采用看似复杂低效的vector+map实现
 *
 * Todo: 把整体变成模板, 学习BGL里面设计property的办法
 * * * * * * * * * * * */
class Neuron{
public:
    const std::vector<const std::string> paras_name_list;
    virtual void Refresh(std::initializer_list<double> input) = 0;
    inline double Measure() const { return spike; }
    virtual ~Neuron() = default;
    std::string PrintUsage() const;
    std::string PrintParaList() const;
    inline double ValueOfParameter(const std::string& para_name){ return parameters[para_name]; }
protected:
    Neuron() = default;
    Neuron(const std::vector<const std::string>& subclass_para_list, std::initializer_list<double> para_inputs);
    double spike;
    std::map<std::string, double> parameters;
};

class Aihara final : public Neuron{
public:
    Aihara(std::initializer_list<double> paras_und_init);
    Aihara(const Aihara& src);
    Aihara& operator=(const Aihara& src);
    Aihara& operator* (double dev);
    void Refresh(std::initializer_list<double> input) override; // refresh neuron status
    inline double MeasureIn() const { return stat_inter; }      // measure the "internal status"
private:
    Aihara() = default;
    double eta_i;        // feedback from other neurons
    double zeta_i;       // 不应性
    double stat_inter;   // internal status
};
/*!
 * @brief Introduction to Hudgkin-Huxley equations
 * mainly refrenced from Hodgkin et.al 1952 in The Journal of Physiology
 * positive current direction is defined from Outside( Extracellular ) to Inside( Intracellular )
 * positive membrane potential direction is defined according to the type of ions:
 * 1. for sodium ions, positive means intracellular potential is higher than extracellular
 * 2. for potassium	ions, positive means extracellular potential is higher than intracellular
 * 3. for 'leakage current' made up by chloride and other ions, it is the same as potassium situation
 *
 * most of the parameter was fixed, only the
 * g_Na and g_K vary with time and membrane potential, decided by
 * g_Na = g_Na_coff * m^3 * h
 * g_K  = g_K_coff  * n^4
 * and time-dependence of gating subunits in the ionic channel:
 *
 * E_L between -49.387 ~ -54.4 mV?
 */
class Hodgkin_huxley final : public Neuron{
public:
    Hodgkin_huxley(std::initializer_list<double> paras_und_init);
    void Refresh(std::initializer_list<double> input) override;
private:
    Hodgkin_huxley() = default;
    // physiological constants
    static const double C_m;        //!< 膜电容(muF/cm^2)
    static const double g_Na_coff;  //!< Na电导(mS/cm^2)
    static const double g_K_coff;   //!< K电导(mS/cm^2)
    static const double g_L;        //!< 漏电导(mS/cm^2)
    static const double E_rest;	    //!< 静息电位(mV) another possible value is -65
    static const double rho_Na;		//!< Na离子浓度(mu*m^-2) 未使用
    static const double rho_K;		//!< Ka离子浓度(mu*m^-2) 未使用

    double E_Na;// = 26 * log(ion_Na_e / ion_Na_i)	mV
    double E_K;	// = 26 * log(ion_K_e / ion_K_i)	mV

    double dV; 	// = E_m - E_rest
    double V_Na, V_K, V_L;
    double m, h, n;
    /* gating variables of ion channels m, h and n
     * (the same dimention of membrane potential).
     */
    inline double alpha_m(double value){ return 0.1*(25.0 + value) / (exp((25.0 + value) / 10.0) - 1); }    // transition rates alpha_m
    inline double alpha_h(double value){ return 0.07*exp(value / 20.0); }               // transition rates alpha_h
    inline double alpha_n(double value){ return 0.01*(value + 10.0) / (exp((value + 10.0) / 10.0) - 1); }   // transition rates alpha_n
    inline double beta_m(double value){ return 4 * exp(value / 18.0); }                 // transition rates beta_m
    inline double beta_h(double value){ return 1 / (exp((30.0 + value) / 10.0) + 1); }  // transition rates beta_h
    inline double beta_n(double value){ return 0.125*exp(value / 80.0); }               // transition rates beta_n
    // voltage - dependent transition rates, which regarded as different
};

class Fitzhugh_nagumo final : public Neuron{
public:
    Fitzhugh_nagumo(std::initializer_list<double> paras_und_init);
    void Refresh(std::initializer_list<double> input) override;
private:
    Fitzhugh_nagumo() = default;
    double w;
};

class Fitzhugh_nagumo2 final : public Neuron{
public:
    Fitzhugh_nagumo2(std::initializer_list<double> paras_und_init);
    void Refresh(std::initializer_list<double> input) override;
private:
    Fitzhugh_nagumo2() = default;
    const double epsilon;
    double w;
};

}   // namespace ntoolkit

#endif //NTOOLKIT_NEURON_H
