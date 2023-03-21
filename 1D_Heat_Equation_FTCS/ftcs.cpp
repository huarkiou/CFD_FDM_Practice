#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <vector>

/* Forward Time Central Space
* ʱ������ǰ��֣��ռ������Ĳ��
*/

int main(void) {
    // �ռ���ɢ
    double x_l = -1.0;
    double x_r = 1.0;
    double dx = 0.025;
    long long nx = static_cast<long long>((x_r - x_l) / dx);

    // ʱ����ɢ
    double dt = 0.0025;
    double t = 1.0;
    long long nt = static_cast<long long>(t / dt);

    double alpha = 1.0 / (M_PI * M_PI);

    // ��ʼ��
    std::vector<double> x(nx + 1), un_cur(nx + 1), un_prev(nx + 1), u_e(nx + 1), error(nx + 1);
    for (long long i = 0; i < nx + 1; ++i) {
        x[i] = x_l + dx * i;
        un_prev[i] = -sin(M_PI * x[i]);
        u_e[i] = -pow(M_E, -t) * sin(M_PI * x[i]);
    }
    un_cur[0] = 0; un_cur[nx] = 0;

    // ������ļ�
    FILE* fp = fopen("./output.txt", "w");
    if (fp == NULL) {
        perror("open output.txt: ");
        exit(1);
    }

    // ����
    double beta = alpha * dt / (dx * dx);

    for (long long k = 0; k < nt + 1; ++k) {
        // ��ָ�ʽ
        for (long long i = 1; i < nx; ++i) {
            un_cur[i] = un_prev[i] + beta * (un_prev[i + 1] - 2 * un_prev[i] + un_prev[i - 1]);
        }
        // �߽�����
        un_cur[0] = 0; un_cur[nx] = 0;
        un_prev = un_cur;
        // �������
        fprintf(fp, "%.6lf\t", k * dt);
        for (int j = 0; j < nx + 1; ++j) {
            fprintf(fp, "%.6lf\t", un_cur[j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    // ����
    double rms_error = 0, max_error = 0;
    for (long long i = 0; i < nx + 1; ++i)
    {
        rms_error = pow(un_cur[i] - u_e[i], 2);
        if (max_error < abs(un_cur[i] - u_e[i])) {
            max_error = abs(un_cur[i] - u_e[i]);
        }
    }
    rms_error = sqrt(rms_error / double(nx + 1));

    printf("rms error = %.8lf, max error = %.8lf\n", rms_error, max_error);
    return 0;
}
