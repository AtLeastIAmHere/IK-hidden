// #include<bits/stdc++.h>
#include<stdio.h>
#include<iostream>
#include<vector>
#include<algorithm>
#include<string>
#include<cstring>
#include<cmath>
#include<set>
#include<fstream>
#include<sstream>
#include<ctime>
#include<map>

#ifdef _WIN64
#include <stdio.h>
#include <io.h>
#include <stdlib.h>
#elif __linux__
#include <sys/types.h>
#include <dirent.h>
#endif

#define EPS (1e-6)

using namespace std;


double stringToDouble(string str) {
    istringstream iss(str);
    double x;
    if (iss >> x)
        return x;
    return 0.0;
} 

string doubleToString(long long d) {
    ostringstream os;
    if (os << d)
        return os.str();
    return "invalid conversion";
}

// long to binary string
string longToBit(long long x){
    string s = "";
    while(x != 0){
        if(x & 1LL) s = '1' + s;
        else s = '0' + s;
        x >>= 1LL;
    }
    return s;
}
 
// #将path中的csv文件读取到strArray中
void readFile(string path, vector<vector<string> > &strArray){
    ifstream inFile(path, ios::in);
    string lineStr;
    
    int flag = 1;
    while(getline(inFile, lineStr)){
        //打印整行字符串
        // cout << lineStr << endl;
        // 存储为二维结构
        stringstream ss(lineStr);
        string str;
        vector<string> lineArray;
        // 按照逗号分隔
        while(getline(ss, str, ','))
            lineArray.push_back(str);
        if(flag){
            flag = 0;
            continue;//跳过表头
        }
        strArray.push_back(lineArray);
    }
}

// 将data以csv的形式写到path中
void writeFile(string &path, vector<vector<int> > &data){
    //静态创建二维数组保存到csv文件
     
    ofstream out(path);
    for(auto &row: data){
        for(int i = 0; i < row.size(); ++i){
            out << row[i];
            if(i != row.size() - 1)
                out << ',';
        }
        out << '\n';
    }
    out.close();
}

// 将data以arff的形式写到path中
void writeArffFile(string &path, vector<vector<int> > &data, string &filename){
    //静态创建二维数组保存到csv文件
     
    // ofstream out(rootPath + resultPath + path10[pi] + ".arff");
    ofstream out(path);

    out << "@RELATION " << filename << endl;
    for(int i = 0; i < data[0].size() - 1; ++i){
        out << "@ATTRIBUTE attribute_" << i << " REAL" << endl;
    }
    out << "@ATTRIBUTE attribute_" << data[0].size() - 1 << " {0,1}" << endl;
    out << endl << "@DATA" << endl;

    for(auto &row: data){
        for(int i = 0; i < row.size(); ++i){
            out << row[i];
            if(i != row.size() - 1)
                out << ',';
        }
        out << '\n';
    }
    out.close();
}

// randomly generate a double number in [0, 1)
double random01(){
    return (rand() % 10000) * 1.0 / 10000;
}


// 检查NDB中是否存在与s0相同的串，如果存在，则说明不匹配
int isMatch(string &s0, vector<set<int> > &NDB){
    for(auto &str: NDB){
        int allBitSame = 1;
        for(auto &x: str){
            // 对于NDB中的每一个串，只要存在一位及以上与s0不同即可
            if(x > 0 && s0[x - 1] == '0') allBitSame = 0;
            if(x < 0 && s0[-x - 1] == '1') allBitSame = 0;
        }

        // NDB中存在与s0一模一样的串，不匹配
        if(allBitSame == 1) return 0; 
    }

    return 1;
}


void getCsvFilesName(string path, vector<string>& files){
    // win和linux读取文件方式略微存在差异
    #ifdef _WIN64
        intptr_t hFile = 0;
        struct _finddata_t fileinfo;
        string p;
        if ((hFile = _findfirst(p.assign(path).append("//*").c_str(), &fileinfo)) != -1){
            do{
                if ((fileinfo.attrib &  _A_SUBDIR)){
                    if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
                        getCsvFilesName(p.assign(path).append("//").append(fileinfo.name), files);
                }
                else{
                    files.push_back(fileinfo.name);
                }
            } while (_findnext(hFile, &fileinfo) == 0);
            _findclose(hFile);
        }

        // 把csv后缀删除
        for(auto &csvName: files){
            csvName.erase(csvName.size() - 4, 4);
        }

    #elif __linux__

        DIR *pDir;
        struct dirent* ptr;
        if(!(pDir = opendir(path.c_str()))){
            cout<<"Folder doesn't Exist!"<<endl;
            return;
        }
        while((ptr = readdir(pDir))!=0) {
            if (strcmp(ptr->d_name, ".") != 0 && strcmp(ptr->d_name, "..") != 0){
                files.push_back(ptr->d_name);
            }
        }
        closedir(pDir);
        for(auto &csvName: files){
            csvName.erase(csvName.size() - 4, 4);
            // cout << csvName << endl;
        }

    #else
        cout << "This code only can run on Windows and Linux!" << endl;
    #endif
}


// csv中存放纯净的数据，无表头；csvName用于返回csv文件名，不含后缀
void dataLoader(string &path, vector<vector<vector<string> > > &data, vector<string> &csvName){

    getCsvFilesName(path, csvName);

    #ifdef _WIN64
        for(auto &csv: csvName){
            vector<vector<string> > csvData;
            string csvPath = path;
            csvPath.append("//").append(csv).append(".csv");

            ifstream inFile(csvPath, ios::in);
            string lineStr;

            while(getline(inFile, lineStr)){
                // 存储为二维结构
                // if(csv == "mysql" || csv == "netbsd") cout << "can get here!" << endl;
                stringstream ss(lineStr);
                string str;
                vector<string> line;
                // 按照逗号分隔
                while(getline(ss, str, ','))
                    line.push_back(str);
                csvData.push_back(line);
            }
            // cout << csv << ":  " << csvData.size() << endl;

            data.push_back(csvData);
        }

    #elif __linux__
        for(auto &csv: csvName){
            vector<vector<string> > csvData;
            string csvPath = path;
            csvPath.append("//").append(csv).append(".csv");

            ifstream inFile(csvPath, ios::in);

            if (!inFile.is_open()) {
                std::cerr << "Error opening file: " << csvPath << std::endl;
                return;
            }

            string lineStr;

            // while(getline(inFile, lineStr)){
            //     vector<string> row;
            //     string cell;
            //     for(char c: lineStr){
            //         if(c == ','){
            //             row.push_back(cell);
            //             cell.clear();
            //         }
            //         else{
            //             cell += c;
            //         }
            //     }
            //     row.push_back(cell);
            //     csvData.push_back(row);
            // }
            while (getline(inFile, lineStr)) {
                // 使用 stringstream 将每行拆分为逗号分隔的值
                istringstream ss(lineStr);
                vector<string> row;

                // 逐个读取每个逗号分隔的值
                string value;
                while (getline(ss, value, ',')) {
                    row.push_back(value);
                }
                // 根据ASCALL码检查行末是否存在回车
                if(static_cast<int>(row[row.size() - 1][row[row.size() - 1].size() - 1]) == 13){
                    // cerr << "Here is a enter in the end of the row!" << endl;
                    row[row.size() - 1] = row[row.size() - 1].substr(0, 1);
                }
                // 将当前行的值添加到 csvData 中
                csvData.push_back(row);
            }

            // while(getline(inFile, lineStr)){
            //     stringstream ss(lineStr);
            //     string str;
            //     vector<string> line;
            //     // 按照逗号分隔
            //     while(getline(ss, str, ','))
            //         line.push_back(str);
            //     csvData.push_back(line);
            // }
            // cout << csv << ":  " << csvData.size() << endl;

            // auto v = csvData[0];
            // int kn = v.size();
            // cout << kn << endl;
            // cout << v[0].size() << " " << v[kn - 2].size() << " " << v[kn - 1].size() << endl;
            // for(char c: v[kn - 1]){
            //     cout << static_cast<int>(c) << " " << endl;
            // }
            // exit(0);

            data.push_back(csvData);
        }
    #else
        cout << "ERROR in function: dataLoader, for the code run on a unknown OS." << endl;
    #endif
}

class FKParameter{
public:
    string s;           // 隐藏串

    int K;              // FK-hidden算法的K值
    vector<double> p;   // 核心K-hidden步骤中的参数，下标0开始
    
    int T;              // 特征数量
    vector<double> f;   // 步骤FK-hidden中的F参数，下标1开始

    int L;              // 特征值在二进制表示下的最大长度
    vector<double> q;   // 步骤QK-hidden中的Q参数，下标1开始

    int m;              // 字符串总长，m = T * L
    int r;              // FK-hidden算法生成的NDB数量阈值参数，算法中生成 N = m * r 个NDB

    double balance = 0; // 范围 [0,1), 平衡处理参数，正值表示过采样，负值表示欠采样
                        // 过采样表示将原始数据中标签少的实例随机复制，直至（标签少的实例数目/总实例数目) >= balance
                        // 欠采样表示将原始数据中标签多的实例随机删除，直至（标签少的实例数目/总实例数目) >= balance

    void init(int _K, int _T, int _L, int _m, int _r){
        this->K = _K;
        this->T = _T;
        this->L = _L;
        this->m = _m;
        this->r = _r;
    }

    void set_balance(double _balance){
        this->balance = _balance;
    }

    void set_p(vector<double> &_p){
        vector <double>().swap(this->p); // 初始化vector
        this->p.resize(_p.size());
        copy(_p.begin(), _p.end(), this->p.begin());
    }

    void set_q(vector<double> &_q){
        vector <double>().swap(this->q);
        this->q.resize(_q.size());
        copy(_q.begin(), _q.end(), this->q.begin());
    }

    void set_f(vector<double> &_f){
        vector <double>().swap(this->f);
        this->f.resize(_f.size());
        copy(_f.begin(), _f.end(), this->f.begin());
    }
};

void getInfoGain(vector<vector<string> > &data_, vector<double> &infoGain, FKParameter &parameter){
    // 每个特征都划分成10个区间，然后根据这些十个区间计算信息增益
    int n = data_.size();
    map<int, double> feature_max, feature_min;
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < parameter.T; ++j){
            feature_max[j] = max(feature_max[j], stringToDouble(data_[i][j]));
            feature_min[j] = min(feature_min[j], stringToDouble(data_[i][j]));
        }
    }
    // 最大最小值没问题

    map<int, map<int, int> > num0, num1; // 记录第i个特征第j个区间有多少个0和1
    map<int, map<int, int> > count_interval;// 记录第i个特征第j个区间有多少实例
    vector<int> label;

    for(int i = 0; i < n; ++i){
        for(int j = 0; j < parameter.T; ++j){
            int x = 10 * (stringToDouble(data_[i][j]) - feature_min[j]) / (feature_max[j] - feature_min[j] + EPS);
            
            // test
            // if(x >= 8)
            // cout << data_[i][j] << "\t" << stringToDouble(data_[i][j]) << "\t" << feature_min[j] << "\t" << feature_max[j] << "\t" << x << endl;

            if(data_[i][parameter.T] == "0") num0[j][x]++;
            else num1[j][x]++;

            count_interval[i][x]++;
        }
    }

    // for(int i = 0; i < parameter.T; ++i){
    //     cout << "feature = " << i << ": ";
    //     for(int j = 0; j < 10; ++j){
    //         cout << count_interval[i][j] << '\t';
    //     }
    //     cout << endl;
    // }
    // exit(0);

    for(int i = 0; i < n; ++i){
        label.push_back(int(stringToDouble(data_[i][parameter.T])));
        // cout << " i = " << i << " label = " << label[i] << "       " << data_[i][parameter.T] << endl;
    }

    // 熵值计算

    // 总熵值E_all
    int num0_all = 0, num1_all = 0;
    for(int i = 0; i < n; ++i){
        if(label[i] == 0) num0_all++;
        else num1_all++;
    }

    double E_all = -1.0 * num0_all / n * log2(1.0 * num0_all / n) - 1.0 * num1_all / n * log2(1.0 * num1_all / n);

    for(int i = 1; i <= parameter.T; ++i) infoGain[i] = E_all;

    // cout << "num1 = " << num1_all << "  num0 =" << num0_all << "   n = " << n << "    n0 / n = " << num0_all / n << "    n1 / n" << num1_all / n << "   log2(n0/n) = " << log2(num0_all / n) << endl;
    // cout << "E_all = " << E_all << endl;exit(0);

    // 计算每个特征的熵值
    for(int i = 0; i < parameter.T; ++i){

        // 计算单个特征中，不同区间的熵值
        for(int interval = 0; interval < 10; ++interval){
            double E_interval = 0;

            if(count_interval[i][interval] == 0) continue;
            if(num0[i][interval] == 0 || num1[i][interval] == 0) continue;

            E_interval -= 1.0 * num0[i][interval] / count_interval[i][interval] * log2(1.0 * num0[i][interval] / count_interval[i][interval]);
            E_interval -= 1.0 * num1[i][interval] / count_interval[i][interval] * log2(1.0 * num1[i][interval] / count_interval[i][interval]);

            infoGain[i + 1] -= 1.0 * count_interval[i][interval] / n * E_interval;
        }
    }    

    // for(int i = 0; i <= parameter.T; ++i){
    //     cout << "i = " << i << ": infoGain = " << infoGain[i] << endl;
    // }
    //         cout << "here" << endl;exit(0);
}

void initial(vector<vector<string> > &data, FKParameter &parameter){


    // T = data[0].size() - 1, minus 1 cause there is a label "bug"
    int K, T, L, r;
    K = 3, T = 82, L = 27, r = 15;
    int m = T * L;

    parameter.init(K, T, L, m, r);


    vector<double> p(parameter.K, 0);
    vector<double> q(parameter.L + 1, 0);

    p[0] = 0.752, p[1] = 0.226, p[2] = 0.022; // p的下标从0开始

    // q的下标从1开始
    for(int i = 1; i <= parameter.L; ++i){
        q[i] = 1.0 / (parameter.L + parameter.L / 2);
    }
    for(int i = 1; i <= parameter.L / 2; ++i){
        q[parameter.L + 1 - i] += 1.0 / (parameter.L + parameter.L / 2);
    }

    //计算不同属性的信息增益，下标从1开始。根据信息增益的占比设置
    vector<double> infoGain(parameter.T + 1, 0);    
    getInfoGain(data, infoGain, parameter);

    // 大于均值的分配更多权重
    double sumInfoGain = 0;
    for(auto &v: infoGain)
        sumInfoGain += v;
    double avg = sumInfoGain / parameter.T;
    for(int i = 1; i <= parameter.T; ++i){
        if(infoGain[i] > avg)
            infoGain[i] = 2;
        else
            infoGain[i] = 1;
        // infoGain[i] /= sumInfoGain;
    }

    // 归一化
    double allWeight = 0;
    for(int i = 1; i <= parameter.T; ++i){
        allWeight += infoGain[i];
    }
    for(int i = 1; i <= parameter.T; ++i){
        infoGain[i] /= allWeight;
    }

    parameter.set_p(p);
    parameter.set_q(q);
    // 使用信息增益的比例初始化FK-hidden的特征选择参数
    parameter.set_f(infoGain);

    srand(unsigned(time(NULL)));
}


// randomly select an attribute t by A[]
int randomAttri(FKParameter &parameter){
    vector<double> preF(parameter.T + 1, 0);
    for(int i = 1; i <= parameter.T; ++i){
        preF[i] = preF[i - 1] + parameter.f[i];
    }
    double rnd = random01();
    for(int i = 1; i <= parameter.T; ++i){
        if(preF[i - 1] <= rnd && rnd < preF[i]){
            return i;
        }
    }

    // no to here
    return 0;
}

// randomly select a bit of the attribure selected by Q[]
int randomBit(FKParameter &parameter){
    vector<double> preQ(parameter.L + 1, 0);
    for(int i = 1; i <= parameter.L; ++i){
        preQ[i] = preQ[i - 1] + parameter.q[i];
    }
    double rnd = random01();
    for(int i = 1; i <= parameter.L; ++i){
        if(preQ[i - 1] <= rnd && rnd < preQ[i]){
            return i;
        }
    }

    // no to here
    return 0;
}

// get the NDB of string s
void getNDB(vector<set<int> > &NDB, FKParameter &parameter){
    int N = parameter.m * parameter.r;
    vector<double> P(parameter.K + 1, 0);
    for(int i = 1; i <= parameter.K; ++i){
        P[i] = P[i - 1] + parameter.p[i - 1];
    }

    while(NDB.size() < N){
        set<int> record;

        double rnd = random01();
        int rndi = 0;
        for(int i = 1; i <= parameter.K; ++i){
            if(P[i - 1] <= rnd && rnd < P[i]){
                rndi = i;
                break;
            }
        }
        // cout << "rndi = " << rndi << endl;
        // cout << "K = " << parameter.K << endl;
        // cout << "rnd = " << rnd << endl;

        // cout << "rnd  = " << rnd << endl;
        // cout << "rndi = " << rndi << endl;

        // generate rndi bits of record which is different from s
        for(int j = 1; j <= rndi; ++j){
            // randomly select an attribute t by A[]
            int attribute_i = randomAttri(parameter);

            // randomly select a bit of the attribure selected by Q[]
            int bit_i = randomBit(parameter);

            // make the selected bit be different from s
            int pos = (attribute_i - 1) * parameter.L + bit_i;
            if(parameter.s[pos - 1] == '1') pos *= -1; // different with s
            // cout << "Diff:   pos = " << pos << "  s[poi] = " << s[abs(pos)] << "   s.size() = " << s.size() << endl;
            record.insert(pos);
        }

        // generate K-rndi bits of record which is same with s for the rest bits with same probability
        while(record.size() < parameter.K){
            int pos = rand() % parameter.m + 1;
            if(record.count(pos) + record.count(-pos) == 0){
                // add rndPos
                if(parameter.s[pos - 1] == '0') pos *= -1;
                record.insert(pos);
            }
        }

        NDB.push_back(record);
        // cout << "record.size() = " << record.size() << endl;

    }// while()
}

void generateTestData(vector<vector<string> > &sourceData, FKParameter &parameter, string &pathTestArff, string &pathTrainCsv, string filename){

    // 1. 归一化数据
    // 2. 取8位小数
    // 3. 8位小数转换为二进制串
    // 4. 二进制串生成NDB，然后生成统计数据
    // 5. 对类不平衡进行处理，进行过采样


    // 1. 数据预处理，归一化数据
    int n1 = sourceData.size();
    int n2 = parameter.T;

    vector<double> minNum, maxNum;
    vector<vector<double> > data;
    vector<vector<string> > data_01;
    vector<int> Y; // class label

    // 记录label
    for(int i = 0; i < n1; ++i){
        Y.push_back(sourceData[i][n2] == "0" ? 0 : 1);
    }

    // 记录每个属性的最大最小值
    for(int i = 0; i < n1; ++i){
        for(int j = 0; j < n2; ++j){
            if(i == 0){
                minNum.push_back(stringToDouble(sourceData[i][j]));
                maxNum.push_back(stringToDouble(sourceData[i][j]));
            }
            else{
                minNum[j] = min(minNum[j], stringToDouble(sourceData[i][j]));
                maxNum[j] = max(maxNum[j], stringToDouble(sourceData[i][j]));
            }
        }
    }

    // 归一化处理
    for(int i = 0; i < n1; ++i){
        vector<double> dataTemp;
        for(int j = 0; j < n2; ++j){
            if(fabs(maxNum[j] - minNum[j]) < 0.000001)
                dataTemp.push_back(0);
            else
                dataTemp.push_back((stringToDouble(sourceData[i][j]) - minNum[j]) / (maxNum[j] - minNum[j])); 
            // skarbonka文件中的noc列，全0，所以此处分母为零
        }
        data.push_back(dataTemp);
    }

    // 使用前8位有效小数
    for(int i = 0; i < n1; ++i){
        for(int j = 0; j < n2; ++j){
            data[i][j] = data[i][j] * 1e8;
        }
    }

    // 转换为二进制
    for(int i = 0; i < n1; ++i){
        vector<string> data_01_temp;
        for(int j = 0; j < n2; ++j){
            string bitString = longToBit((long long)(data[i][j]));
            // 补充数据长度
            //十进制8位数，转换为二进制的位数最多是L(L=27)位
            while(bitString.size() < parameter.L){
                bitString = "0" + bitString;
            }            
            data_01_temp.push_back(bitString);
        }
        data_01.push_back(data_01_temp);
    }

    // NDB 01 字符统计
    vector<vector<int> > extraction;
    for(int i = 0; i < n1; ++i){

        // generate NDB
        string s0 = "";
        for(int j = 0; j < n2; ++j)
            s0 = s0 + data_01[i][j];
        parameter.s = s0;
        // cout << s << endl << s.size() << endl;
        vector<set<int> > NDB;
        getNDB(NDB, parameter);
        if(isMatch(s0, NDB) == 0){
            cout << "!!!!!NDB is not matched!!!!!" << endl;
            exit(0);
        }

        // feature extraction 
        vector<int> num(2 * parameter.m, 0);
        for(auto &line: NDB){
            for(auto &x: line){
                int pos  = x > 0 ? x - 1 : -x - 1 + parameter.m;
                // cout << "x = " << x << "  pos = " << pos << endl;
                num[pos]++;                
            }
        }
        extraction.push_back(num);
    // cout << "i can get here!" << endl; exit(0);
    }

    // add class label
    for(int i = 0; i < n1; ++i){
        extraction[i].push_back(Y[i]);
        if(Y[i] != 0 && Y[i] != 1){
            cout << "n1 = " << n1 << endl;
            cout << "line:" << i << endl;
            exit(0);
        }
    }

    // save data
    writeFile(pathTrainCsv, extraction);

    // // write to arff file
    writeArffFile(pathTestArff, extraction, filename);
}


void generateTrainData(vector<vector<string> > &sourceData, FKParameter &parameter, string &pathTrainArff, string &pathTrainCsv, string filename){

    // 1. 归一化数据
    // 2. 取8位小数
    // 3. 8位小数转换为二进制串
    // 4. 二进制串生成NDB，然后生成统计数据
    // 5. 对类不平衡进行处理，进行过采样


    // 1. 数据预处理
    int n1 = sourceData.size();
    int n2 = n1 > 0 ? sourceData[0].size() - 1 : 0;
    
    vector<double> minNum, maxNum;
    vector<vector<double> > data;
    vector<vector<string> > data_01;
    vector<int> Y; // class label

    // 记录label
    for(int i = 0; i < n1; ++i){
        // cout << sourceData[i][0] << "  " << sourceData[i][n2] << endl;
        Y.push_back(sourceData[i][n2] == "0" ? 0 : 1);
    }

    // 数据平衡
    int num_defective = 0;
    int num_nondefective = 0;
    vector<int> less_id_set;
    vector<int> more_id_set;
    for(int i = 0; i < n1; ++i){
        if(Y[i] == 0) num_nondefective++;
        else num_defective++;
    }

    // 记录数目少和数量更多的标签是什么
    int label = num_defective < num_nondefective ? 1 : 0;
    int num_min = min(num_defective, num_nondefective);
    int num_max = max(num_defective, num_nondefective);
    for(int i = 0; i < n1; ++i){
        if(Y[i] == label) 
            less_id_set.push_back(i);
        else 
            more_id_set.push_back(i);
    }

    // 过采样，从已有数据集中类别少的数据随机选取，进行复制
    if(parameter.balance > 0){
        // int num_add = parameter.balance * num_max - (1.0 - parameter.balance) * num_min;
        int num_add = (parameter.balance * (num_max + num_min) - num_min) / (1.0 - parameter.balance);
        // cout << "num_add = " << num_add << endl;
        // cout << less_id_set.size() << "   " << more_id_set.size() << endl;
        // exit(0);
        for(int i = 1; i <= num_add; ++i){
            int randId = rand() % less_id_set.size();
            sourceData.push_back(sourceData[less_id_set[randId]]);
        }

    }

    // 欠采样，将已有数据集中类别更多的数据，随机删除，知道满足要求
    if(parameter.balance < 0){
        int num_delete = 1.0 * num_max - num_min * (1.0 - parameter.balance) / (parameter.balance);
        if(num_delete < -10000){
            cout << "a = " << num_min << "   b = " << num_max << "   m = " << parameter.balance << endl;
        }
        // cout << "num_delete = " << num_delete << endl;
        map<int, int> id_delete; // 标记需要删除的id
        for(int i = 1; i <= num_delete; ++i){
            while(1){
                int randId = rand() % more_id_set.size();
                if(id_delete[more_id_set[randId]] == 0){
                    id_delete[more_id_set[randId]] = 1;
                    break;
                } 
            }
        }

        // 将标记的id对应的元素删除，因为从前删除会影响后面的id，所以需要从后往前删除
        for(int i = sourceData.size() - 1; i >= 0; --i){
            if(id_delete[i] == 1){
                sourceData.erase(sourceData.begin() + i);
            }
        }
    }

    // cout << "  Balance = " << parameter.balance << endl;
    // cout << "Add " << int(int(sourceData.size()) - n1) << " new records\n";

    // 将过拟合新增数据的标签记录下来
    for(int i = n1; i < sourceData.size(); ++i){
        Y.push_back(sourceData[i][n2] == "0" ? 0 : 1);
    }
    n1 = sourceData.size();
        
    // cout << "n1 = " << n1 << "   " << sourceData.size() << endl;
    // cout << "n2 = " << n2 << "   " << sourceData[0].size() << endl;

    // 记录每个属性的最大最小值
    for(int i = 0; i < n1; ++i){
        for(int j = 0; j < n2; ++j){
            if(i == 0){
                minNum.push_back(stringToDouble(sourceData[i][j]));
                maxNum.push_back(stringToDouble(sourceData[i][j]));
            }
            else{
                minNum[j] = min(minNum[j], stringToDouble(sourceData[i][j]));
                maxNum[j] = max(maxNum[j], stringToDouble(sourceData[i][j]));
            }
        }
    }

    // 归一化处理
    for(int i = 0; i < n1; ++i){
        vector<double> dataTemp;
        for(int j = 0; j < n2; ++j){
            if(fabs(maxNum[j] - minNum[j]) < 0.000001)
                dataTemp.push_back(0);
            else
                dataTemp.push_back((stringToDouble(sourceData[i][j]) - minNum[j]) / (maxNum[j] - minNum[j])); 
            // skarbonka文件中的noc列，全0，所以此处分母为零
        }
        data.push_back(dataTemp);
    }

    // 使用前8位有效小数
    for(int i = 0; i < n1; ++i){
        for(int j = 0; j < n2; ++j){
            data[i][j] = data[i][j] * 1e8;
        }
    }

    // 转换为二进制
    for(int i = 0; i < n1; ++i){
        vector<string> data_01_temp;
        for(int j = 0; j < n2; ++j){
            string bitString = longToBit((long long)(data[i][j]));
            // 补充数据长度
            //十进制8位数，转换为二进制的位数最多是L(L=27)位
            while(bitString.size() < parameter.L){
                bitString = "0" + bitString;
            }            
            data_01_temp.push_back(bitString);
        }
        data_01.push_back(data_01_temp);
    }


    // NDB 01 字符统计
    vector<vector<int> > extraction;
    for(int i = 0; i < n1; ++i){

        // generate NDB
        string s0 = "";
        for(int j = 0; j < n2; ++j)
            s0 = s0 + data_01[i][j];
        parameter.s = s0;
        // cout << s << endl << s.size() << endl;
        vector<set<int> > NDB;
        getNDB(NDB, parameter);
        if(isMatch(s0, NDB) == 0){
            cout << "!!!!!NDB is not matched!!!!!" << endl;
            exit(0);
        }

        // feature extraction 
        vector<int> num(2 * parameter.m, 0);
        for(auto &line: NDB){
            for(auto &x: line){
                int pos  = x > 0 ? x - 1 : -x - 1 + parameter.m;
                // cout << "x = " << x << "  pos = " << pos << endl;
                num[pos]++;                
            }
        }
        extraction.push_back(num);
    }

    // add class label
    for(int i = 0; i < n1; ++i){
        extraction[i].push_back(Y[i]);
        if(Y[i] != 0 && Y[i] != 1){
            cout << Y[i] << endl;
            cout << "line:" << i << endl;
            exit(0);
        }
    }

    // save data as csv
    writeFile(pathTrainCsv, extraction);

    // save data as arff
    writeArffFile(pathTrainArff, extraction, filename);
}


// sourceData中包含多个表格，将每个表格生成一个训练集，一个测试集
// 对于表格a.csv，生成训练集train_a.arff，train_a.csv，测试集test_a.arff，test_a.csv。
void makeSingleData(vector<vector<vector<string> > > &sourceData, vector<string> &csvName, string &pathResult, double balance){
    for(int i_train; i_train < csvName.size(); ++i_train){

        cout << "[MAKING SINGLE DATA]------Running train data: " << csvName[i_train] << endl;

        FKParameter parameter;
        initial(sourceData[i_train], parameter);// 根据原始数据，计算信息增益的熵值，根据熵值对FK-hidden算法的参数进行初始化
        parameter.set_balance(balance);    // 类不平衡参数

        for(int i_test = 0; i_test < csvName.size(); ++i_test){
            if(i_test == i_train) continue;

            cout << "                 Running test data: " << csvName[i_test] << endl;

            // 文件名如 train_file1TOfile2，表示使用file1作为训练集，且是只用于file2测试的训练集
            string fileName = csvName[i_train] + "TO" + csvName[i_test];
            string pathTrainArff = pathResult + "//train_" + fileName + ".arff";
            string pathTrainCsv = pathResult + "//train_" + fileName + ".csv";
            generateTrainData(sourceData[i_train], parameter, pathTrainArff, pathTrainCsv, "train_" + fileName);

            // 生成测试集，利用训练集的InfoGain参数进行生成
            string pathTestArff = pathResult + "//test_" + fileName + ".arff";
            string pathTestCsv = pathResult + "//test_" + fileName + ".csv";
            generateTestData(sourceData[i_test], parameter, pathTestArff, pathTestCsv, "test_" + fileName);
            cout << "                 finish test data: " << csvName[i_test] << endl;
        }
        cout << "[MAKING SINGLE DATA]------Finish train data: " << csvName[i_train] << endl << endl;
    }
}

void makeCrossData(vector<vector<vector<string> > > &sourceData, vector<string> &csvName, string &pathResult, double balance){

    for(int i = 0; i < csvName.size(); ++i){
        cout << "[MAKING CROSS DATA]------Running data: " << csvName[i] << endl;

        FKParameter parameter;
        vector<vector<string> > dataTrain;
        for(int j = 0; j < csvName.size(); ++j){
            if(i == j) continue;
            for(auto &v: sourceData[j]){
                dataTrain.push_back(v);
            }
        }

        initial(dataTrain, parameter);// 根据原始数据，计算信息增益的熵值，根据熵值对FK-hidden算法的参数进行初始化
        parameter.set_balance(balance);    // 类不平衡参数

        string pathTrainArff = pathResult + "//train_" + csvName[i] + ".arff";
        string pathTrainCsv = pathResult + "//train_" + csvName[i] + ".csv";
        generateTrainData(dataTrain, parameter, pathTrainArff, pathTrainCsv, "train_" + csvName[i]);

        // 生成测试集，利用训练集的InfoGain参数进行生成
        string pathTestArff = pathResult + "//test_" + csvName[i] + ".arff";
        string pathTestCsv = pathResult + "//test_" + csvName[i] + ".csv";
        generateTestData(sourceData[i], parameter, pathTestArff, pathTestCsv, "test_" + csvName[i]);

        cout << "[MAKING CROSS DATA]------Finish csv file: " << csvName[i] << endl << endl;
    }
}


int main(){
    string pathSource = "data//source"; // 存放csv文件的目录
    vector<string> csvName;
    vector<vector<vector<string> > > sourceData; // 将三个表格存储到sourceData中
    dataLoader(pathSource, sourceData, csvName);

    // string pathResultSingle = "data//resultSingle";
    // string pathResultCross = "data//resultCross";

    // 重复10次实验
    int repeatTimes = 10;
    for(int idx_repeat = 0; idx_repeat < repeatTimes; ++idx_repeat){

        // 类不平衡处理参数，正值表示过采样，负值表示欠采样
        for(int i = 0; i <= 0; i += 5){
            string pathResultSingle = string("data//result//").append(doubleToString(idx_repeat)).append("//").append(doubleToString(i)).append("//resultSingle");
            string pathResultCross = string("data//result//").append(doubleToString(idx_repeat)).append("//").append(doubleToString(i)).append("//resultCross");
            // cout << pathResultSingle << endl; exit(0);
            double balance = 1.0 * i / 100;

            cout << endl << "[RUNNING BALANCE]------" << balance << endl;
            // makeSingleData(sourceData, csvName, pathResultSingle, balance); // 单个数据集生成测试集，生成一个训练集
            makeCrossData(sourceData, csvName, pathResultCross, balance);   // 选择其中一个数据集作为测试集，剩下合并为一个训练集
            cout << "[FINISH BALANCE]------" << balance << endl << endl;
        } 
    }
    return 0;
}