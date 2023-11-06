
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void calc_sin_repeat(float x, float out[100]) {
    out[0] = sin(x);
    out[1] = sin(x);
    out[2] = sin(x);
    out[3] = sin(x);
    out[4] = sin(x);
    out[5] = sin(x);
    out[6] = sin(x);
    out[7] = sin(x);
    out[8] = sin(x);
    out[9] = sin(x);
    out[10] = sin(x);
    out[11] = sin(x);
    out[12] = sin(x);
    out[13] = sin(x);
    out[14] = sin(x);
    out[15] = sin(x);
    out[16] = sin(x);
    out[17] = sin(x);
    out[18] = sin(x);
    out[19] = sin(x);
    out[20] = sin(x);
    out[21] = sin(x);
    out[22] = sin(x);
    out[23] = sin(x);
    out[24] = sin(x);
    out[25] = sin(x);
    out[26] = sin(x);
    out[27] = sin(x);
    out[28] = sin(x);
    out[29] = sin(x);
    out[30] = sin(x);
    out[31] = sin(x);
    out[32] = sin(x);
    out[33] = sin(x);
    out[34] = sin(x);
    out[35] = sin(x);
    out[36] = sin(x);
    out[37] = sin(x);
    out[38] = sin(x);
    out[39] = sin(x);
    out[40] = sin(x);
    out[41] = sin(x);
    out[42] = sin(x);
    out[43] = sin(x);
    out[44] = sin(x);
    out[45] = sin(x);
    out[46] = sin(x);
    out[47] = sin(x);
    out[48] = sin(x);
    out[49] = sin(x);
    out[50] = sin(x);
    out[51] = sin(x);
    out[52] = sin(x);
    out[53] = sin(x);
    out[54] = sin(x);
    out[55] = sin(x);
    out[56] = sin(x);
    out[57] = sin(x);
    out[58] = sin(x);
    out[59] = sin(x);
    out[60] = sin(x);
    out[61] = sin(x);
    out[62] = sin(x);
    out[63] = sin(x);
    out[64] = sin(x);
    out[65] = sin(x);
    out[66] = sin(x);
    out[67] = sin(x);
    out[68] = sin(x);
    out[69] = sin(x);
    out[70] = sin(x);
    out[71] = sin(x);
    out[72] = sin(x);
    out[73] = sin(x);
    out[74] = sin(x);
    out[75] = sin(x);
    out[76] = sin(x);
    out[77] = sin(x);
    out[78] = sin(x);
    out[79] = sin(x);
    out[80] = sin(x);
    out[81] = sin(x);
    out[82] = sin(x);
    out[83] = sin(x);
    out[84] = sin(x);
    out[85] = sin(x);
    out[86] = sin(x);
    out[87] = sin(x);
    out[88] = sin(x);
    out[89] = sin(x);
    out[90] = sin(x);
    out[91] = sin(x);
    out[92] = sin(x);
    out[93] = sin(x);
    out[94] = sin(x);
    out[95] = sin(x);
    out[96] = sin(x);
    out[97] = sin(x);
    out[98] = sin(x);
    out[99] = sin(x);
}
               
void calc_sin_repeat_for_loop(float x, float out[100]) {
    static int i;
    for(i=0;i<100;i++) {
        out[i] = sin(x);
    }
}
               
void calc_sin_no_repeat(float x, float out[100]) {
    int i;
    float sin_x = sin(x);
    for(i=0;i<100;i++) {
        out[i] = sin_x;
    }
}
               
               
void test_swap() {
    float a[] = {1,2,5};
    float b[] = {4,8,2};
    
    // create a pointer that points to b array
    float* b_ptr = b;
    // create a pointer that points to a array
    float* a_ptr = a;

    printf("a[0] = %f\n",a_ptr[0]);
    printf("a[1] = %f\n",a_ptr[1]);
    printf("a[2] = %f\n",a_ptr[2]);
    printf("b[0] = %f\n",b_ptr[0]);
    printf("b[1] = %f\n",b_ptr[1]);
    printf("b[2] = %f\n",b_ptr[2]);

    float* tmp = b_ptr;
    b_ptr = a_ptr;
    a_ptr = tmp;

    printf("test_swap2\n");
    printf("a[0] = %f\n",a_ptr[0]);
    printf("a[1] = %f\n",a_ptr[1]);
    printf("a[2] = %f\n",a_ptr[2]);
    printf("b[0] = %f\n",b_ptr[0]);
    printf("b[1] = %f\n",b_ptr[1]);
    printf("b[2] = %f\n",b_ptr[2]);
}
