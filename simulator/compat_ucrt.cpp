// El libraylib.a incluido se compiló con un MinGW basado en msvcrt, donde stat() se
// resuelve a stat64i32. Los MinGW basados en UCRT solo exportan _stat64i32, así que
// aquí se da el nombre que falta. CMake solo añade este archivo si hace falta.
extern "C" {
struct _stat64i32;
int _stat64i32(const char* path, struct _stat64i32* buffer);
int stat64i32(const char* path, struct _stat64i32* buffer) { return _stat64i32(path, buffer); }
}
