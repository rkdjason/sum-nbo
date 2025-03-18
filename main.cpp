#include <stdio.h>
#include <stdint.h>

uint32_t my_ntohl(uint32_t n){
        return (n >> 24) | (n & 0xFF0000) >> 8 | (n & 0xFF00) << 8 | (n << 24);
}

int main(int argc, char *argv[]){
        if (argc < 2) {
                printf("Syntax : %s <file1> [<file2>...]\n", argv[0]);
                return 1;
        }

        uint32_t sum = 0;
        for(int i = 1; i < argc; i++){
                FILE *fp = NULL;
                long int fsize;
                uint32_t n;

                if ((fp = fopen(argv[i], "rb")) == NULL) {
                        printf("Error : file %s is not opened\n", argv[i]);
                        return 1;
                }

                fseek(fp, 0, SEEK_END);
                fsize = ftell(fp);
                if (fsize != 4){
                        printf("Error : size of file %s is not 4 bytes\n", argv[i]);
                        fclose(fp);
                        return 1;
                }
                fseek(fp, 0, SEEK_SET);

                if (fread(&n, sizeof(n), 1, fp) != 1){
                        printf("Error : can not read file %s\n", argv[i]);
                        fclose(fp);
                        return 1;
                }
                fclose(fp);

                n = my_ntohl(n);
                sum += n;
                printf("%u(0x%04x)", n, n);
                if (i < argc - 1){
                        printf(" + ");
                }
        }

        printf(" = %u(0x%04x)\n", sum, sum);

        return 0;
}
