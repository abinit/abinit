#define ABITRACE(msg) write(0,'("=====DEBUG===== ",I4," in file ",A , "in function", A, "message: ",A )') __LINE__,__FILE__, __FUNC__, msg
#define ABITRACE2(msg) write(msg, '("file: "A, " line:" , I4)', __FILE__, __LINE__)

