#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

extern "C" {

int makedir(const char * path){
    struct stat st = {0};

    if (stat(path,&st) == -1) {
        return mkdir(path,0700);
    } else {
        return 0;
    }
}

}
