#define _USE_MATH_DEFINES
#define PXL_SIZE 512
#define BLK_SIZE 8
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

unsigned char raw_arr[PXL_SIZE][PXL_SIZE];
unsigned char block[64][64][BLK_SIZE][BLK_SIZE];
float block_flt[64][64][BLK_SIZE][BLK_SIZE];
char scanned[64][64][64];
unsigned char out_arr[PXL_SIZE][PXL_SIZE];
const int BLK_NUM = PXL_SIZE / BLK_SIZE;
char quant[8][8] = {
    {16,11,10,16,24,40,51,61},
    {12,12,14,19,26,58,60,55},
    {14,13,16,24,40,57,69,56},
    {14,17,22,29,51,87,80,62},
    {18,22,37,56,68,109,103,77},
    {24,35,55,64,81,104,113,92},
    {49,64,78,87,103,121,120,101},
    {72,92,95,98,112,100,103,99} };

void Block() {
    for (int i = 0; i < BLK_NUM; i++) {
        for (int j = 0; j < BLK_NUM; j++) {
            for (int k = 0; k < BLK_SIZE; k++) {
                for (int l = 0; l < BLK_SIZE; l++)
                    block[i][j][k][l] = raw_arr[((i * BLK_SIZE) + k)][((j * BLK_SIZE) + l)];
            }
        }
    }
    printf("Blocked\n");
}
void DCT() {
    float cv, cu;
    int i, j, v, u, y, x = 0;
    for (i = 0; i < BLK_NUM; i++) {
        for (j = 0; j < BLK_NUM; j++) {
            for (v = 0; v < BLK_SIZE; v++) {
                for (u = 0; u < BLK_SIZE; u++) {
                    if (v == 0)
                        cv = 1 / sqrt(2);
                    else cv = 1;
                    if (u == 0)
                        cu = 1 / sqrt(2);
                    else cu = 1;

                    for (y = 0; y < BLK_SIZE; y++) {
                        for (x = 0; x < BLK_SIZE; x++)
                            block_flt[i][j][v][u] += 0.25 * cv * cu * (float)block[i][j][y][x] * 
                                cos((((2 * y + 1) * v * M_PI)) / 16) * cos((((2 * x + 1) * u * M_PI)) / 16);
                    }
                }
            }
        }
    }
    printf("DCT\n");
}
void Quant() {
    int i, j, v, u = 0;
    for (i = 0; i < BLK_NUM; i++) {
        for (j = 0; j < BLK_NUM; j++) {
            for (v = 0; v < BLK_SIZE; v++) {
                for (u = 0; u < BLK_SIZE; u++)
                    block_flt[i][j][v][u] = roundf(block_flt[i][j][v][u] / quant[v][u]);
            }
        }
    }
    printf("Quantization\n");
}
void ZigZag() {
    int i, j, y, x, m = 0;
    for (i = 0; i < BLK_NUM; i++) {
        for (j = 0; j < BLK_NUM; j++) {
            int index = 0;
            for (m = 0; m < BLK_SIZE; m++) {
                if (m % 2 == 0) {
                    y = m;
                    x = 0;
                    while (y >= 0) {
                        scanned[i][j][index] = (char)block_flt[i][j][y][x];
                        y--;
                        x++;
                        index++;
                    }
                }
                else {
                    y = 0;
                    x = m;
                    while (x >= 0) {
                        scanned[i][j][index] = (char)block_flt[i][j][y][x];
                        y++;
                        x--;
                        index++;
                    }
                }
            }
            for (m = 1; m < BLK_SIZE; m++) {
                if (m % 2 != 0) {
                    y = BLK_SIZE - 1;
                    x = m;
                    while (y >= m) {
                        scanned[i][j][index] = (char)block_flt[i][j][y][x];
                        y--;
                        x++;
                        index++;
                    }
                }
                else {
                    y = m;
                    x = BLK_SIZE - 1;
                    while (x >= m) {
                        scanned[i][j][index] = (char)block_flt[i][j][y][x];
                        y++;
                        x--;
                        index++;
                    }
                }
            }
        }
    }
    printf("Zig-Zag Scanning\n");
}
void BZigZag() {
    int i, j, y, x, m = 0;
    for (i = 0; i < BLK_NUM; i++) {
        for (j = 0; j < BLK_NUM; j++) {
            int index = 0;
            for (m = 0; m < BLK_SIZE; m++) {
                if (m % 2 == 0) {
                    y = m;
                    x = 0;
                    while (y >= 0) {
                        block_flt[i][j][y][x] = (float)scanned[i][j][index];
                        y--;
                        x++;
                        index++;
                    }
                }
                else {
                    y = 0;
                    x = m;
                    while (x >= 0) {
                        block_flt[i][j][y][x] = (float)scanned[i][j][index];
                        y++;
                        x--;
                        index++;
                    }
                }
            }
            for (m = 1; m < BLK_SIZE; m++) {
                if (m % 2 != 0) {
                    y = BLK_SIZE - 1;
                    x = m;
                    while (y >= m) {
                        block_flt[i][j][y][x] = (float)scanned[i][j][index];
                        y--;
                        x++;
                        index++;
                    }
                }
                else {
                    y = m;
                    x = BLK_SIZE - 1;
                    while (x >= m) {
                        block_flt[i][j][y][x] = (float)scanned[i][j][index];
                        y++;
                        x--;
                        index++;
                    }
                }
            }
        }
    }
    printf("\nReconstruction\n");
}
void IQuant() {
    int i, j, v, u = 0;
    for (i = 0; i < BLK_NUM; i++) {
        for (j = 0; j < BLK_NUM; j++) {
            for (v = 0; v < BLK_SIZE; v++) {
                for (u = 0; u < BLK_SIZE; u++)
                    block_flt[i][j][v][u] *= quant[v][u];
            }
        }
    }
    printf("Inv. Quantization\n");
}
void IDCT() {
    float cv, cu;
    int i, j, v, u, y, x = 0;
    for (i = 0; i < BLK_NUM; i++) {
        for (j = 0; j < BLK_NUM; j++) {
            for (y = 0; y < BLK_SIZE; y++) {
                for (x = 0; x < BLK_SIZE; x++) {
                    float out_sum = 0;
                    for (v = 0; v < BLK_SIZE; v++) {
                        for (u = 0; u < BLK_SIZE; u++) {
                            if (v == 0)
                                cv = 1 / sqrt(2);
                            else cv = 1;
                            if (u == 0)
                                cu = 1 / sqrt(2);
                            else cu = 1;
                            out_sum += 0.25 * cv * cu * block_flt[i][j][v][u] * cos((((2 * y + 1) * v * M_PI)) / 16) * cos((((2 * x + 1) * u * M_PI)) / 16);
                            block[i][j][y][x] = (unsigned char)roundf(out_sum);
                        }
                    }
                }
            }
        }
    }
    printf("IDCT\n");
}
void Unblock() {
    for (int i = 0; i < BLK_NUM; i++) {
        for (int j = 0; j < BLK_NUM; j++) {
            for (int k = 0; k < BLK_SIZE; k++) {
                for (int l = 0; l < BLK_SIZE; l++)
                    out_arr[((i * BLK_SIZE) + k)][((j * BLK_SIZE) + l)] = block[i][j][k][l];
            }
        }
    }
    printf("Unblocked\n");
}
float MSE(unsigned char* raw[], unsigned char* out[]) {
    float err_rate = 0;
    for (int i = 0; i < PXL_SIZE; i++) {
        for (int j = 0; j < PXL_SIZE; j++)
            err_rate += pow(raw_arr[i][j] - out_arr[i][j], 2);
    }
    err_rate /= (PXL_SIZE * PXL_SIZE);

    return err_rate;
}
void JPEG() {
    Block();
    DCT();
    Quant();
    ZigZag();
    BZigZag();
    IQuant();
    IDCT();
    Unblock();
    printf("\nMSE Error Rate : %.2f%\n", MSE(raw_arr, out_arr));
}

int main() {
    FILE* raw = NULL;
    FILE* out = NULL;
    fopen_s(&raw, "C:\\Users\\deep\\Desktop\\JPEG\\lena.raw", "r");
    fopen_s(&out, "C:\\Users\\deep\\Desktop\\JPEG\\lena_out.raw", "w");

    if (raw == NULL || out == NULL) {
        printf("Failed to fetch/write the file!");
        return 1;
    }
    fread(raw_arr, 1, PXL_SIZE * PXL_SIZE, raw);
    JPEG();
    fwrite(out_arr, 1, PXL_SIZE * PXL_SIZE, out);

    fclose(raw);
    fclose(out);

    return 0;
}