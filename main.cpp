/*!
 * @file main.cpp
 * @brief 这是一个神经网络工具箱
 * @details 没什么细节
 * @mainpage Neuron Network Toolkit
 * @author Tooko Madobe
 * @email madobet@outlook.com
 * @version 0.1
 * @date 2019-09-05
 * @license GPLv3
 */

#include <unistd.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <random>
#include <aio.h>

#include "lib.h"
#include "neuron.h"
#include "network.h"
#include "model.h"

#ifdef LINUX
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#elif WIN
#include <io.h>
#include <direct.h>
#endif

//#define FIXED_STAT_DETECT // check whether all the output is fixed, if it is, the terminate this simulation and abort
using namespace std;
class globalArgs_t {
public:
    globalArgs_t():neuron_names({ "Aihara","Hodgkin-Huxley","FitzHugh-Nagumo","FitzHugh-Nagumo2" }),
                   net_names({ "brain-net" })
                   {
        net_size = 0;
        neuron_name.clear();//
        net_name.clear();  // clear model name string
        para.clear();       // clear the parameter vector
        dt = 0.1;           // simulation interval, units: ms, ref: spike usually last for 1 ms
        time = 1200.0;      // total time, units: ms
        abort_time = time/2;

        is_rand = true;      // should there be is_rand initialization for each neuron
        init_file.clear();  // specify the initialization file name to empty string
        init_val = 0.0;     // make the initial value to be zero

        w_file.clear();//
        w_folder="weight";
        w_thd = 0.2;  //
        cp_str = 1.0;        //

        sti_func.clear();   //

        cfg_file.clear();//
        out_file.clear();   //
        verbosity = 0;
    }
    ~globalArgs_t(){}
    int ModID(const char* neuron_name, const char* net_name)
    /* 根据给出的模型名称给出相应的编号，
     * 编号是由构造函数根据模型名称的字符串在neuron_mods_name、net_mods_name中的顺序自动产生的，
     * 存储在mods_n变量中，
     * 规则是：神经元模型编号为n，网络模型编号为m，则模型编号为1000*m+n
     * 调用时应当自己决定查询的是神经元还是网络模型
     * Q：那么为什么不直接用枚举？
     * A：因为我希望编号的同时，也能查询模型的名字，枚举只能通过名字得到编号，反过来不行
     *    另外因为对全部组合都进行编号，所以用枚举全部要手动编写反而不方便
     * Q：那么为什么不直接用vector
     * A：也是相同的道理
     * Q：multimap呢？
     * A：还没试过
     * 最好用哈希实现！
     * */
    {
        int id = -1;
        for(int i=0; i < neuron_names.size(); i++){
            if(neuron_names.at(i) == neuron_name){
                id = i;
                break;
            }
        }
        for(int i=0; i < net_names.size(); i++){
            if(net_names.at(i) == net_name){
                id += i*1000;
                break;
            }
        }
        return id;
    }
    vector<const char*> ModName(int id) const
    // ModID的逆操作
    {
        int neuron_id = id%1000;
        int net_id = id/1000;
        return vector<const char*>{ neuron_names.at(neuron_id), net_names.at(net_id) };
    }

    const vector<const char*> neuron_names;
    const vector<const char*> net_names;
    int net_size;
    string neuron_name;     /* -n option, neuron model name     */
    string net_name;        /* -N option, net model name    */
    vector<double> para;    /* -P option, model parameters      */
    double dt;              /* -d option, ms                    */
    double time;            /* -t option, ms                    */
    double abort_time;      /* -a */
    // model option

    bool is_rand;
    string init_file;       /* -I option    */
    double init_val;        /* -i option    */
    // initialization option
    // the initial order is is_rand > initial file > initial value

    string w_file;          /* -W option, weight file name */
    string w_folder;        /* -F option, */
    double w_thd;
    double cp_str;          /* -c option, coupled strength */
    // weight option

    string sti_func;        /* -f option, extarnal stimulate function, can only accept liner or sin, I_stim uA? */
    // stimulation option
    // stimulate order sti file > offset+last+str+pat

    string cfg_file;         /* -C option */
    string out_file;         /* -O option, output file name */
    int verbosity;           /* -v option,  */
    // program option
}glb_args;
// 全局参数

bool ProcessOptions(int agc, char* agv[]);
bool CheckOptions();
void DisplayUsage();
void ConvertDocument();
int main(int argc, char *argv[])
{
    if(!ProcessOptions(argc, argv)){
        return -1;
    };
    if(!CheckOptions()){
        return -2; // Check option
    }

//---------------------------------------------------------------------------------------
    // typedef <typename T><vector<T>> matrix<T>; 类似这样的要如何表达？
    vector< vector< vector<double> > > weight_list; weight_list.clear();
    if(glb_args.w_folder.empty()){
        try{
            weight_list.push_back(ntoolkit::io::ReadMatFromFile<double>(glb_args.w_file));
        }catch(const char* msg){
            cerr << msg << endl;
        }
    }else{
#ifdef LINUX
        DIR *dir;
        struct dirent *dir_ptr;
        string folder_path("./");
        folder_path += glb_args.w_folder;    // 目录下不可包含子目录！
        dir = opendir(folder_path.c_str());
        while((dir_ptr = readdir(dir)) != nullptr){
            string file_path; file_path.clear();
            if(!strncmp(dir_ptr->d_name,".",128)||!strncmp(dir_ptr->d_name,"..",128))
                continue;
            file_path = folder_path + "/" + dir_ptr->d_name;
            try{
                weight_list.push_back(ntoolkit::io::ReadMatFromFile(file_path, glb_args.w_thd));
            }catch(const char* msg){
                cerr << msg << endl;
            }
            //
        }
        closedir(dir);
#endif
    }
    // Import weight;

    vector<double> initial; initial.clear();
    if(glb_args.verbosity > 0) cout << "Initial Voltage:\n";
    if(glb_args.is_rand){
        static default_random_engine e;
        e.seed(time(nullptr));
        if(glb_args.neuron_name == glb_args.neuron_names.at(0)){
            uniform_real_distribution<double> u(0,1); // Generate between 0 ~ 1
            for(int i = 0; i < glb_args.net_size; i++){
                initial.push_back(u(e));
                if(glb_args.verbosity > 0)  cout << initial.back() << "   ";
            }
        }else if(glb_args.neuron_name == glb_args.neuron_names.at(1)){
            uniform_real_distribution<double> u(-70,35); // Generpaate between -70 ~ 35
            for(int i = 0; i < glb_args.net_size; i++){
                initial.push_back(u(e));
                if(glb_args.verbosity > 0)  cout << initial.back() << "   ";
            }
        }
    } else if(!glb_args.init_file.empty()){
        ifstream ifs(glb_args.init_file);
        if(!ifs.is_open()) {
            cerr << "ERROR:initial file not found" << endl;
            ifs.close();
            return -3;
        }
        for(int i = 0; i < glb_args.net_size; i++){
            double temp;
            ifs >> temp;
            initial.push_back(temp);
            if(glb_args.verbosity > 0)  cout << temp << "   ";
        }
        ifs.close();
    } else {
        for(int i = 0; i < glb_args.net_size; i++){
            initial.push_back(glb_args.init_val);
            if(glb_args.verbosity > 0)  cout << glb_args.init_val << "   ";
        }
    }
    if(glb_args.verbosity > 0) cout << endl;
    // Import initial vector

	const double abnormal_per = 0.1; // Percentage of abnormal value can be tolerated.
    auto abnormal_show = int(glb_args.time/glb_args.dt*abnormal_per);
    int abnormal_cnt = abnormal_show;
#ifdef FIXED_STAT_DETECT
    const double max_spike_inter = 125; // ms
    int samedata_show = int(max_spike_inter/glb_args.dt);
    int samedata_cnt = samedata_show;
    // same data can not continously over samedata_show, ALpha wave min = 8hz, so the maxmium between two spike should no more than 125 ms
    double dif_per = 0.05;  // how much percent of the difference is allowed, below which is considered to be the same.
#endif
    if(glb_args.verbosity > 0){
        cout \
		<< "# Experiment starts in " << glb_args.time << " ms with " \
		<< glb_args.dt << " ms/step using " \
		<< glb_args.neuron_name <<  glb_args.net_name <<"model.\n"
		<< "# Result will be export to " << glb_args.out_file \
		<< endl;
    }

    vector<vector<vector<double>>> all_samples_spike;
    vector<vector<double>> all_sample_act;

    vector<double> all_sample_act_ave;
    vector<double> all_sample_act_min;
    vector<double> all_sample_act_max;
    for(int i=0; i<glb_args.net_size; i++){
        all_sample_act_ave.push_back(0);
        all_sample_act_min.push_back(0);
        all_sample_act_max.push_back(0);
    }
    const double act_unit = glb_args.dt/(glb_args.time-glb_args.abort_time);
    const int sample_n = weight_list.size();

    for( auto &it : weight_list ){
        ntoolkit::Neuron_net_base *model = nullptr;
        if(glb_args.net_name == glb_args.neuron_names.at(0)){
            if(glb_args.neuron_name == glb_args.neuron_names.at(0)){
                model = new ntoolkit::Neuron_net<ntoolkit::Aihara, ntoolkit::Brain_net>(initial, it, {
                    glb_args.para.at(0),glb_args.para.at(1),glb_args.para.at(2),glb_args.para.at(3)
                    });
            }else if(glb_args.neuron_name == glb_args.neuron_names.at(1)){
                model = new ntoolkit::Neuron_net<ntoolkit::Hodgkin_huxley, ntoolkit::Brain_net>(initial, it, {
                        glb_args.para.at(0),glb_args.para.at(1),glb_args.para.at(2),glb_args.para.at(3),glb_args.para.at(4),glb_args.dt
                });
            }else if(glb_args.neuron_name == glb_args.neuron_names.at(3)){
                model = new ntoolkit::Neuron_net<ntoolkit::Fitzhugh_nagumo2, ntoolkit::Brain_net>(initial, it, {
                        glb_args.para.at(0),glb_args.para.at(1),glb_args.para.at(2),glb_args.para.at(3),glb_args.dt
                });
            }else{
                return -5;
            }
        }

        model->GenNetReport("aaa.txt");
        vector<vector<double>> single_sample_spike;
        vector<vector<double>> single_sample_sti;
        vector<double> res_t0;
        vector<double> single_sample_act;
        //
        single_sample_spike.clear();
        single_sample_sti.clear();
        single_sample_act.clear();
        stringstream ss;    // converter
        ss.clear();
        string time_str;    // store time
        double time = 0.0;

        for(int i=0; i<glb_args.net_size; i++){
            res_t0.push_back(0);
            single_sample_act.push_back(0);
        }
        while ( time < glb_args.time )
        {
            vector<double> ext; ext.clear();
            if(!glb_args.sti_func.empty()){
                for( int i=0; i < glb_args.net_size; i++ ){

                    ss.clear();
                    time_str.clear();
                    ss << time;
                    ss >> time_str;
                    // convert double to string

                    double temp=ntoolkit::math::CalFunc(glb_args.sti_func,"x",time_str);
                    if(glb_args.verbosity > 0) cout << temp << endl;
                    ext.push_back(temp);
                }
                single_sample_sti.push_back(ext);
            }else{ // 若sti为空，则直接用方
                double local_time = time-(int)time/4*4;
                for( int i=0; i < glb_args.net_size; i++ ){
                    if(local_time >= 0 && local_time < 0.5){
                        ext.push_back(0.0);
                    }else if(local_time >= 0.5 && local_time < 4){
                        ext.push_back(0.75);
                    }
                }
                single_sample_sti.push_back(ext);
            }// Generate stimulation

            vector<double> res_t1;
            res_t1 = model->Refresh(ext);
            if(time > glb_args.abort_time) {
                for (int i=0; i<glb_args.net_size; i++) {
                    if(!isnormal(res_t1.at(i))) {
                        //abnormal_cnt--;
                    }
                    if(abs(res_t1.at(i)-res_t0.at(i)) >= 0.5){
                        single_sample_act[i] += act_unit;
                        all_sample_act_ave[i] += act_unit/sample_n;
                    }
                }
                if (abnormal_cnt < 0) {
                    cerr << "ERROR:invalid neuron spike more than "<< abnormal_show*glb_args.dt <<" ms" << endl;
                    //delete model;
                    return -4;
                }
#ifdef FIXED_STAT_DETECT
                if(!single_sample_spike.empty()) {
            vector<double> re = single_sample_spike.back();
            bool sameflag = true;
            for (int i = 0; i < res_t1.size(); i++) {
                if ( fabs((res_t1.at(i) - re.at(i))/re.at(i)) > dif_per ) {
                    sameflag = false;
                    break;
                }
            }
            if(sameflag) {
                samedata_cnt--;
            }
            else
                samedata_cnt = samedata_show;
        }
        if(samedata_cnt < 0){
            cerr << "ERROR:similar value appears last more than "<< samedata_show*glb_args.dt <<" ms, differnece less than " << dif_per*100 << "\%" << endl;
            return -5;
        }// 3rd check out
#endif
                single_sample_spike.push_back(res_t1);
            }
            for(int i=0; i<glb_args.net_size; i++){
                if(single_sample_act.at(i) > all_sample_act_max.at(i)){
                    all_sample_act_max[i] = single_sample_act.at(i);
                }
                if(single_sample_act.at(i) < all_sample_act_min.at(i)){
                    all_sample_act_min[i] = single_sample_act.at(i);
                }
            }
            res_t0 = res_t1;
            time = time + glb_args.dt;
        }// Simulation
        all_samples_spike.push_back(single_sample_spike);
        all_sample_act.push_back(single_sample_act);

        delete model;
    }

    /*
     * 输出结果
     * */
    string save_file_name;
    string para_str;
    ostringstream convert;convert.clear();
    for(const auto &it : glb_args.para){
        convert << it;
        convert << '_';
    }
    para_str = convert.str();
    for(int i=0; i<weight_list.size(); i++){
        convert.clear();
        convert << "sample_";
        convert << i;
        save_file_name = glb_args.neuron_name + '_' + glb_args.net_name + '_' + para_str + convert.str() + ".txt";
        save_file_name.clear();
        ntoolkit::io::WriteMat2File(all_samples_spike[i], save_file_name);
    }
    vector<vector<double>> res2;
    save_file_name = glb_args.neuron_name + '_' + glb_args.net_name + para_str + "_sample_all_ave-min-max.txt";
    res2.clear();
    res2.push_back(all_sample_act_ave);
    res2.push_back(all_sample_act_min);
    res2.push_back(all_sample_act_max);
    ntoolkit::io::WriteMat2File(res2, save_file_name);
    save_file_name = glb_args.neuron_name + '_' + glb_args.net_name + para_str + "_sample_all_ave-min-max_left.txt";
    res2.clear();
    ntoolkit::io::WriteMat2File(res2, save_file_name);
    save_file_name = glb_args.neuron_name + '_' + glb_args.net_name + para_str + "_sample_all_ave-min-max_right.txt";
    res2.clear();
    ntoolkit::io::WriteMat2File(res2, save_file_name);
    save_file_name = glb_args.neuron_name + '_' + glb_args.net_name + para_str + "_sample_all.txt";
    ntoolkit::io::WriteMat2File(all_sample_act, save_file_name);

    // interested areas
    int select_area;
    select_area = 1;
    vector<vector<double>> res;
    for( const auto &it: weight_list){
        res.push_back(ExtractOneAreaSpike(it,select_area-1));
    }
    save_file_name = glb_args.neuron_name + '_' + glb_args.net_name + para_str + "_area_1_spike_5000_5500.txt";
    ntoolkit::io::WriteMat2File(res, save_file_name);
    select_area = 2;
    res.clear();
    for( const auto &it: weight_list){
        res.push_back(ExtractOneAreaSpike(it,select_area-1));
    }
    save_file_name = glb_args.neuron_name + '_' + glb_args.net_name + para_str + "_area_2_spike_5000_5500.txt";
    ntoolkit::io::WriteMat2File(res, save_file_name);

    select_area = 25;
    res.clear();
    for( const auto &it: weight_list){
        res.push_back(ExtractOneAreaSpike(it,select_area-1));
    }
    save_file_name = glb_args.neuron_name + '_' + glb_args.net_name + para_str + "_area_25_spike_5000_5500.txt";
    ntoolkit::io::WriteMat2File(res, save_file_name);
    select_area = 26;
    res.clear();
    for( const auto &it: weight_list){
        res.push_back(ExtractOneAreaSpike(it,select_area-1));
    }
    save_file_name = glb_args.neuron_name + '_' + glb_args.net_name + para_str + "_area_26_spike_5000_5500.txt";
    ntoolkit::io::WriteMat2File(res, save_file_name);

    select_area = 79;
    res.clear();
    for( const auto &it: weight_list){
        res.push_back(ExtractOneAreaSpike(it,select_area-1));
    }
    save_file_name = glb_args.neuron_name + '_' + glb_args.net_name + para_str + "_area_79_spike_5000_5500.txt";
    ntoolkit::io::WriteMat2File(res, save_file_name);
    select_area = 80;
    res.clear();
    for( const auto &it: weight_list){
        res.push_back(ExtractOneAreaSpike(it,select_area-1));
    }
    save_file_name = glb_args.neuron_name + '_' + glb_args.net_name + para_str + "_area_80_spike_5000_5500.txt";
    ntoolkit::io::WriteMat2File(res, save_file_name);


    // 权重分布情况
    res.clear();
	return 0;
}

bool ProcessOptions(int agc, char* agv[]){
    const char *optString = "n:N:P:d:t:a:I:i:W:F:c:l:r:f:C:O:vh?";
    stringstream converter;
    int opt = getopt( agc, agv, optString );
    while( opt != -1 ) {
        switch( opt ){
            case 'n':
                glb_args.neuron_name = optarg;
                break;
            case 'N':
                glb_args.net_name = optarg;
                break;
            case 'P':
                converter.clear();
                converter.str(optarg);
                while(true){
                    double temp;
                    converter >> temp;
                    if( converter.fail() ) break;
                    glb_args.para.push_back(temp);
                }
                break;
            case 'd':
                converter.clear();
                converter.str(optarg);
                {
                    double temp;
                    converter >> temp;
                    glb_args.dt = temp;
                }
                break;
            case 't':
                converter.clear();
                converter.str(optarg);
                {
                    double temp;
                    converter >> temp;
                    glb_args.time = temp;
                }
                break;
            case 'a':
                converter.clear();
                converter.str(optarg);
                {
                    double temp;
                    converter >> temp;
                    glb_args.abort_time = temp;
                }
                break;
                // Neuron_network_base
            case 'I':
                glb_args.is_rand = false;
                glb_args.init_file = optarg;
                break;
            case 'i':
                glb_args.is_rand = false;
                converter.clear();
                converter.str(optarg);
                {
                    double temp;
                    converter >> temp;
                    glb_args.init_val = temp;
                }
                break;
                // Initial
            case 'W':
                glb_args.w_file = optarg;
                break;
            case 'F':
                glb_args.w_folder = optarg;
                break;
            case 'c':
                converter.clear();
                converter.str(optarg);
                {
                    double temp;
                    converter >> temp;
                    glb_args.cp_str = temp;
                }
                break;
                // Weight
            case 'f':
                glb_args.sti_func = optarg;
                break;
            case 'C':
                glb_args.cfg_file = optarg;
                break;
            case 'O':
                glb_args.out_file = optarg;
                break;
            case 'v':
                glb_args.verbosity++;
                break;
            case 'h':
            case '?':
                DisplayUsage();
                break;
            default:
                /* You won't actually get here. */
                break;
        }
        opt = getopt( agc, agv, optString );
        converter.clear();
    }
    return true;
}
bool CheckOptions() {
    if ( glb_args.neuron_name.empty() ) {
        cerr << "ERROR:neuron model unspecified" << endl;
        return false;
    }else if( glb_args.net_name.empty() ){
        cerr << "ERROR: net model unspecified" << endl;
        return false;
    }else{
        bool not_found = true;
        for( const string& cp : glb_args.neuron_names ){
            if(glb_args.neuron_name == cp) {
                not_found = false;
                break;
            }
        }
        if(not_found){
            cerr << "ERROR:unknown neuron model" << endl;
            return false;
        }
        not_found = true;
        for( const string& cp : glb_args.net_names ){
            if(glb_args.net_name == cp) {
                not_found = false;
                break;
            }
        }
        if(not_found){
            cerr << "ERROR:unknown net model" << endl;
            return false;
        }
    }
    if( glb_args.time < 0.0 ){
        cerr << "ERROR: time must bigger than zero" << endl;
        return false;
    }
    if( glb_args.abort_time < 0.0 ){
        cerr << "ERROR: abort time must bigger than zero" << endl;
        return false;
    }else if( glb_args.abort_time > glb_args.time ) {
        cerr << "ERROR: abort time bigger than time" << endl;
        return false;
    }
    if( glb_args.w_folder.empty() ) {
        cerr << "INFO: weight folder unspecified, try to use weight file" << endl;
    }else if( glb_args.w_file.empty() ) {
        cerr << "ERROR: no folder no file, weight matrix unloaded, supposed to be single neuron" << endl;
        return false;
    }
    if( glb_args.para.empty() ){
        cout << "WARN: model without parameter" << endl;
    } else if( glb_args.cp_str == 0 ) {
        cout << "WARN: net is uncoupled" << endl;
    }
    if( !glb_args.is_rand && glb_args.init_file.empty() &&glb_args.init_val == 0.0 ){
        cout << "WARN: all neurons initialized with 0.0" << endl;
    }
    if( glb_args.sti_func.empty() ) {
        cout << "INFO: no stimulate function specified, without stimulation" << endl;
    }
    if( !glb_args.cfg_file.empty() ) {
        cout << "INFO: using config file instead, ignore all options" << endl;
    }
    if( glb_args.out_file.empty() ) {
        cout << "INFO: using default out file name" << endl;
        glb_args.out_file = "result.txt";
    }
    return true;
}
void DisplayUsage(){
    cout\
    << "\nmodel options:\n"\
    << "\t-n\tspecify a neuron model\n"\
    << "\t-N\tspecify a net model\n"\
    << "\t-P\tparameter list\n"\
    << "\t-d\tdt\n"\
    << "\t-t\ttime\n"\
    << "\ninitialization options:\n"\
    << "\t-I\timport initial status from a file\n"\
    << "\t-i\tinitial value\n"\
    << "\nWeight options:\n"\
    << "\t-W\timport weight from a file\n"\
    << "\t-F\timport weight from a folder\n"\
    << "\t-c\tcoupled strength\n"\
    << "\nStimulation option:\n"\
    << "\t-f\tfunction\n"\
    << "\nother options:\n"\
    << "\t-C\tconfig file\n"\
    << "\t-O\toutput file name\n"\
    << "\t-v\tverbosity mode"\
    << endl;
    exit( EXIT_FAILURE );
}
// Convert the input files to HTML, governed by glb_args.
void ConvertDocument()
{
    /* ... */
}
/*
int CStr2Int(const char* str){
}

double CStr2Double(const char* str){
}
*/
