#include <iostream>
#include <algorithm>
#include <fstream>
#include <iomanip>
const double Pi = acos(-1.0);//圆周率
const double EARTH_RADIUS = 6378.137;//计算距离所需常数
const double my_eps = 1e-16;//极小值，用来判0
int recur_num;//用来计算递归层数的全局变量

class eNodeB{//基站类
public:
    int id;//基站id
    double lon;//经度
    double lat;//纬度
    double k_dist;//用来求第k小的k_dist
    eNodeB(int idx=0,double lonx=0,double latx=0,double distx=0) {//构造函数
        id = idx;
        lon = lonx;
        lat = latx;
        k_dist = distx; 
    }
    eNodeB(const eNodeB& b) {//拷贝构造函数
        id = b.id;
        lon = b.lon;
        lat = b.lat;
        k_dist = b.k_dist;
    }
    void get(int idx = 0, double lonx = 0, double latx = 0, double distx = 0) {//赋值函数
        id = idx;
        lon = lonx;
        lat = latx;
        k_dist = distx;
    }
    double operator-(const eNodeB& b) {//两个基站类相减返回值为距离
        double radlon = lon * Pi / 180;
        double radlat = lat * Pi / 180;
        double radlatb = b.lat * Pi / 180;
        double radlonb= b.lon * Pi / 180;
        double v = acos(cos(radlat) * cos(radlatb) * cos(radlon - radlonb) + sin(radlat) * sin(radlatb));
        v = v * EARTH_RADIUS;
        v *= 1000;
        if (v >= 0)return v;
        else return -v;
    }
    void operator=(const eNodeB& b) {//赋值运算符重载
        id = b.id;
        lon = b.lon;
        lat = b.lat;
        k_dist = b.k_dist;
    }
};
eNodeB f1, f2;//最近的两个点
eNodeB s1, s2;//次近的两个点
double ans1 = -1;//最近的距离，未赋值时初始化为负数
double ans2 = -1;//次近的距离，未赋值时初始化为负数
/* int Get_eNodeB(eNodeB* &arr)
* 从文件中读入基站数据
* 输入数组头指针的引用
*/
int Get_eNodeB(eNodeB* &arr) {
    std::ifstream enodebTXT;
    enodebTXT.open("enodeb.txt", std::ios::in);
    int n;
    enodebTXT >> n;
    arr = new eNodeB[n + 10];
    int id;
    double lon, lat, dist;
    for (int i = 0; i < n; ++i) {
        enodebTXT >> id >> lon >> lat >> dist;
        arr[i].get(id, lon, lat, dist);
    }
    return n;
}
/* eNodeB tri_Get_K_dist(eNodeB* arr, int left, int right, int k)
* 二分实现线性时间选择
* arr为数组头指针，left和right为区间两端下标，k表示查找的值为当前区间第k小
*/
eNodeB Get_K_dist(eNodeB* arr, int left, int right, int k, int cur) {
    recur_num = cur > recur_num ? cur : recur_num;
    if (left >= right)return arr[left];
    int l = left, h = right;
    int id = left + rand() % (right - left + 1);
    eNodeB temp = arr[id];
    arr[id] = arr[left];
    arr[left] = temp;
    for (; l < h;) {
        while (l < h && arr[h].k_dist >= temp.k_dist) {
            --h;
        }arr[l] = arr[h];
        while (l < h && arr[l].k_dist <= temp.k_dist) {
            ++l;
        }arr[h] = arr[l];
    }arr[l] = temp;
    if (l > k) return Get_K_dist(arr, left, l - 1, k, cur + 1);
    else if (l < k) return Get_K_dist(arr, l + 1, right, k, cur + 1);
    return temp;
}
/* eNodeB Tri_Get_K_dist(eNodeB* arr, int left, int right, int k)
* 三分实现线性时间选择
* arr为数组头指针，left和right为区间两端下标，k表示查找的值为当前区间第k小
*/
eNodeB Tri_Get_K_dist(eNodeB* arr, int left, int right, int k, int cur) {
    recur_num = cur > recur_num ? cur : recur_num;
    //std::cout <<recur_num<<' '<< left << ' ' << right << std::endl;
    if (left >= right)return arr[left];
    int id = left + rand() % (right - left + 1);
    eNodeB temp = arr[id], z;
    arr[id] = arr[left];
    arr[left] = temp;
    int l = left, h = right, i = left + 1;
    for (; i <= h;) {
        if(arr[i].k_dist > temp.k_dist) {
            z = arr[h];
            arr[h] = arr[i];
            arr[i]=z;
            --h;
        }
        else if(arr[i].k_dist < temp.k_dist) {
            z = arr[l];
            arr[l] = arr[i];
            arr[i] = z;
            ++l; ++i;
        }
        else {
            ++i;
        }
    }
    if (l > k) return Tri_Get_K_dist(arr, left, l - 1, k, cur + 1);
    else if (h < k) return Tri_Get_K_dist(arr, h + 1, right, k, cur + 1);
    return temp;
}
bool cmplat(eNodeB x, eNodeB y) { //以纬度排序
    if(x.lat!=y.lat)return x.lat<y.lat;
    return x.lon < y.lon;
}
bool cmplon(eNodeB x, eNodeB y) { //以经度排序
    if (x.lon != y.lon)return x.lon < y.lon;
    return x.lat < y.lat;
}
/* inline bool judge_ans(int idx, int idy)
* 判断当前点是否已经在答案里，避免计算两次。
* 参数为两个基站的坐标
*/
inline bool judge_ans(int idx, int idy) { 
    if (idx == f1.id && idy == f2.id)return 1;
    if (idx == f2.id && idy == f1.id)return 1;
    if (idx == s1.id && idy == s2.id)return 1;
    if (idx == s2.id && idy == s1.id)return 1;
    return 0;
}
/* void upd_ans(eNodeB x, eNodeB y)
* 用当前点更新最近和次近平面点对的答案
* arr为数组头指针，left和right为区间两端下标，
* temp为临时处理需要的额外空间
*/
void upd_ans(eNodeB x, eNodeB y) {
    if (judge_ans(x.id, y.id))return;
    double v = x - y;
    if (v < my_eps) return;
    if (ans1<0 || ans1>v) {
        ans2 = ans1; s1 = f1; s2 = f2;
        ans1 = v; f1 = x; f2 = y;
    }
    else if (ans2<0 || ans2>v) {
        ans2 = v; s1 = x; s2 = y;
    }
}
/* void my_merge(eNodeB* arr, int left, int right, eNodeB* temp)
* 递归二分计算平面最近点对
* arr为数组头指针，left和right为区间两端下标，
* temp为临时处理需要的额外空间
*/
void my_merge(eNodeB* arr, int left, int right, eNodeB* temp, int cur) {
    recur_num = cur > recur_num ? cur : recur_num;
    if (right - left <= 3) {
        for (int i = left; i <= right; ++i)
            for (int j = i + 1; j <= right; ++j) upd_ans(arr[i], arr[j]);
        std::sort(arr + left, arr + right + 1, cmplat);
        return;
    }
    int mid = (left + right) / 2;
    double midlon = arr[mid].lon;
    my_merge(arr, left, mid, temp, cur + 1);
    my_merge(arr, mid + 1, right, temp, cur + 1);
    for (int i = left, j = mid + 1, w = left; w <= right; ) {
        if (i <= mid && (j > right || arr[i].lat < arr[j].lat)) {
            temp[w++] = arr[i++];
        }
        else {
            temp[w++] = arr[j++];
        }
    }
    for (int i = left; i <= right; ++i) arr[i] = temp[i];
    int tn = 0;
    for (int i = left; i <= right; ++i) {
        if (abs(arr[i].lon - midlon) < ans1) {
            for (int j = tn - 1; j >= 0 && arr[i].lat - temp[j].lat < ans1; --j)
                upd_ans(arr[i], temp[j]);
            temp[tn++] = arr[i];
        }
    }
}
/* void GetPairPoints(eNodeB* arr,int siz)
* 计算平面最近点对
* arr为数组头指针，siz为数组大小
*/
void GetPairPoints(eNodeB* arr,int siz) {
    std::sort(arr,arr+siz,cmplon);
    eNodeB* temp = new eNodeB[siz + 10];
    my_merge(arr, 0, siz - 1, temp, 1);
}
int main() {
    std::fstream answerTXT;
    answerTXT.open("answer.txt", std::ios::out);
    srand((unsigned int)time(NULL));
    eNodeB* arr=NULL;
    int siz;
    siz=Get_eNodeB(arr);
    std::cout << "读入数据完毕" << std::endl;
    recur_num = 0;
    int test[4] = {1,5,50,1033};
    for (int i = 0; i < 4; ++i) {
        int k = test[i];
        recur_num = 0;
        eNodeB ans = Get_K_dist(arr, 0, siz - 1, k - 1, 1);
        answerTXT << "第" << k << "小值的基站id为" << ans.id << std::endl;
        answerTXT << "第" << k << "小值为" << std::fixed<< std::setprecision(3) << ans.k_dist << std::endl;
        answerTXT << "二分查找循环层数为" << recur_num << std::endl;
        answerTXT << std::endl;

        recur_num = 0;
        ans = Tri_Get_K_dist(arr, 0, siz - 1, k - 1, 1);
        answerTXT << "第" << k << "小值的基站id为" << ans.id << std::endl;
        answerTXT << "第" << k << "小值为" << std::fixed << std::setprecision(3) << ans.k_dist << std::endl;
        answerTXT << "三分查找循环层数为" << recur_num << std::endl;
        answerTXT << std::endl;
    }

    recur_num = 0;
    GetPairPoints(arr, siz);
    answerTXT << "循环层数为" << recur_num << std::endl;
    answerTXT << "最近两点的id为" << f1.id<<' '<<f2.id << std::endl;
    answerTXT << "最小距离为" << std::fixed << std::setprecision(3) << ans1 << std::endl;
    answerTXT << "次近两点的id为" << s1.id << ' ' << s2.id << std::endl;
    answerTXT << "次小距离为" << std::fixed << std::setprecision(3) << ans2 << std::endl;
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
