#include "VtkopContract.h"
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkNew.h>
#include <vtkNonOverlappingAMR.h>
#include <vtkPoints.h>
#include <vtkTetra.h>
#include <vtkUniformGrid.h>
#include <vtkImageData.h>
#include <vtkVersionMacros.h>
#include <vtkUnstructuredGrid.h>
#include <vtkAMRDualContour.h>
#include <vtkAppendPolyData.h>
#include <vtkCompositeDataIterator.h>
#include <vtkDataObject.h>
#include <vtkDummyController.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiProcessController.h>
#include <vtkPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkCommand.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
namespace {
void check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
class Errors final : public vtkCommand {
public:
    static Errors* New(){return new Errors;}
    bool failed=false;
    void Execute(vtkObject*,unsigned long,void*) override {failed=true;}
};
double exportDC(vtkNonOverlappingAMR* amr,const std::filesystem::path& path){
    struct ControllerScope {
        vtkSmartPointer<vtkMultiProcessController> previous;
        vtkNew<vtkDummyController> serial;
        ControllerScope():previous(vtkMultiProcessController::GetGlobalController()){
            vtkMultiProcessController::SetGlobalController(serial);
        }
        ~ControllerScope(){vtkMultiProcessController::SetGlobalController(previous);}
    } scope;
    vtkNew<Errors> errors;
    vtkNew<vtkAMRDualContour> contour;
    contour->AddObserver(vtkCommand::ErrorEvent,errors);
    contour->SetInputDataObject(amr);
    contour->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_CELLS,"density");
    contour->SetIsoValue(.5);
    contour->SetEnableCapping(true);
    contour->SetEnableMergePoints(true);
    contour->SetController(scope.serial);
    contour->SetEnableMultiProcessCommunication(1);
    const auto dcStart=std::chrono::steady_clock::now();
    contour->Update();
    const double dcSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-dcStart).count();
    check(!errors->failed,"AMRDualContour failed");
    auto output=vtkMultiBlockDataSet::SafeDownCast(contour->GetOutputDataObject(0));
    check(output!=nullptr,"Missing DC output");
    vtkNew<vtkAppendPolyData> append;
    auto iterator=vtkSmartPointer<vtkCompositeDataIterator>::Take(output->NewIterator());
    iterator->SkipEmptyNodesOn();
    int blocks=0;
    for(iterator->InitTraversal();!iterator->IsDoneWithTraversal();iterator->GoToNextItem()){
        auto poly=vtkPolyData::SafeDownCast(iterator->GetCurrentDataObject());
        if(poly&&poly->GetNumberOfPolys()>0){append->AddInputData(poly);++blocks;}
    }
    check(blocks>0,"DC produced an empty surface");
    vtkNew<vtkTriangleFilter> triangles;
    triangles->AddObserver(vtkCommand::ErrorEvent,errors);
    triangles->SetInputConnection(append->GetOutputPort());
    triangles->PassLinesOff();triangles->PassVertsOff();triangles->Update();
    check(!errors->failed,"DC triangulation failed");auto mesh=triangles->GetOutput();
    check(mesh->GetNumberOfPolys()>0,"Empty triangulated surface");
    std::ofstream obj(path);check(bool(obj),"Cannot open OBJ");obj.precision(17);
    for(vtkIdType i=0;i<mesh->GetNumberOfPoints();++i){double p[3];mesh->GetPoint(i,p);check(std::isfinite(p[0])&&std::isfinite(p[1])&&std::isfinite(p[2]),"Nonfinite DC vertex");obj<<"v "<<p[0]<<' '<<p[1]<<' '<<p[2]<<'\n';}
    for(vtkIdType i=0;i<mesh->GetNumberOfCells();++i){auto c=mesh->GetCell(i);check(c->GetNumberOfPoints()==3,"Nontriangle DC face");obj<<"f";for(int j=0;j<3;++j)obj<<' '<<c->GetPointId(j)+1;obj<<'\n';}
    obj.flush();check(bool(obj),"OBJ write failed");
    std::cout<<"OBJ: "<<path.string()<<" ("<<mesh->GetNumberOfPoints()<<" vertices, "<<mesh->GetNumberOfPolys()<<" triangles)\n";
    return dcSeconds;
}
bool inside(const std::string& shape,const double* p){
    const double x=p[0],y=p[1],z=p[2];
    if(shape=="sphere")return x*x+y*y+z*z<.65*.65;
    if(shape=="cone")return z>=-.5 && z<=.5 && std::hypot(x,y)<.25-.5*z;
    if(shape=="torus")return std::hypot(std::hypot(x,y)-.48,z)<.18;
    constexpr double k=6.283185307179586;
    return std::max({std::abs(x),std::abs(y),std::abs(z)})<.78 &&
        std::abs(std::sin(k*x)*std::cos(k*y)+std::sin(k*y)*std::cos(k*z)+std::sin(k*z)*std::cos(k*x))<.35;
}
vtkSmartPointer<vtkUnstructuredGrid> makeSource(const std::string& shape,int n){
    vtkNew<vtkPoints> points;points->SetDataTypeToDouble();
    auto id=[n](int x,int y,int z)->vtkIdType{return (vtkIdType(z)*(n+1)+y)*(n+1)+x;};
    const double h=2./n;
    for(int z=0;z<=n;++z)for(int y=0;y<=n;++y)for(int x=0;x<=n;++x){
        double p[3]={-1+x*h,-1+y*h,-1+z*h};
        if(x>0&&x<n&&y>0&&y<n&&z>0&&z<n)
            for(int d=0;d<3;++d)p[d]+=.06*h*std::sin(1.7*x+2.3*y+3.1*z+1.2*d);
        points->InsertNextPoint(p);
    }
    auto mesh=vtkSmartPointer<vtkUnstructuredGrid>::New();mesh->SetPoints(points);
    vtkNew<vtkDoubleArray> vf;vf->SetName("vf");
    for(int z=0;z<n;++z)for(int y=0;y<n;++y)for(int x=0;x<n;++x){
        std::array<int,3> permutation{0,1,2};
        do{
            int q[3]={x,y,z};vtkIdType ids[4];ids[0]=id(q[0],q[1],q[2]);
            for(int j=0;j<3;++j){++q[permutation[j]];ids[j+1]=id(q[0],q[1],q[2]);}
            double p[4][3],center[3]={0,0,0};
            for(int j=0;j<4;++j){points->GetPoint(ids[j],p[j]);for(int d=0;d<3;++d)center[d]+=.25*p[j][d];}
            double v=vtkTetra::ComputeVolume(p[0],p[1],p[2],p[3]);
            check(std::isfinite(v)&&std::abs(v)>1e-15,"Degenerate source tetrahedron");
            if(v<0)std::swap(ids[1],ids[2]);
            if(inside(shape,center)){mesh->InsertNextCell(VTK_TETRA,4,ids);vf->InsertNextValue(1.);}
        }while(std::next_permutation(permutation.begin(),permutation.end()));
    }
    mesh->GetCellData()->SetScalars(vf);check(mesh->GetNumberOfCells()>0,"Empty source");return mesh;
}
double sourceVolume(vtkUnstructuredGrid* mesh){double sum=0;for(vtkIdType i=0;i<mesh->GetNumberOfCells();++i){double p[4][3];auto c=mesh->GetCell(i);for(int j=0;j<4;++j)c->GetPoints()->GetPoint(j,p[j]);sum+=std::abs(vtkTetra::ComputeVolume(p[0],p[1],p[2],p[3]));}return sum;}
}
int main(int argc,char** argv){try{
    check(argc==2,"Usage: VFTransferShapes <new-output-directory>");
    const std::filesystem::path out=argv[1];check(!std::filesystem::exists(out),"Output directory must not exist");std::filesystem::create_directories(out);
    std::ofstream timings(out/"dc_timings.txt");check(bool(timings),"Cannot open DC timings");
    timings<<"DC wall time: vtkAMRDualContour::Update only, single invocation per case.\n"
           <<"Excludes VF transfer, triangulation and OBJ writing. Not a repeated benchmark.\n"
           <<"source blocks: 64^3; target: 150 cells per source AABB axis; no padding; iso: 0.5\n\n";
    timings<<std::fixed<<std::setprecision(6);
    constexpr int sourceN=64,targetN=150;
    for(const std::string shape:{"sphere","gyroid","torus","cone"}){
        auto source=makeSource(shape,sourceN);const double reference=sourceVolume(source);
        const auto count=source->GetNumberOfCells();const auto pointCount=source->GetNumberOfPoints();
        const auto grid=makeTransferGrid(source,{targetN,targetN,targetN});
        for(const auto method:{TransferMethod::CellCenterSampling,TransferMethod::NodalInterpolation,TransferMethod::VolumeWeighted}){
            const char* name=method==TransferMethod::CellCenterSampling?"center":method==TransferMethod::NodalInterpolation?"nodal":"volume";
            const auto start=std::chrono::steady_clock::now();auto amr=transferVF(source,"vf",grid,method);
            const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            check(amr->GetNumberOfLevels()==1,"Wrong AMR levels");
#if VTK_VERSION_NUMBER >= VTK_VERSION_CHECK(9,6,0)
            auto block=amr->GetDataSetAsImageData(0,0);
#else
            auto block=amr->GetDataSet(0,0);
#endif
            check(block!=nullptr,"Missing AMR block");
            auto rho=block->GetCellData()->GetArray("density");
            check(rho,"Missing density array");check(rho->GetNumberOfTuples()==vtkIdType(grid.cells[0])*grid.cells[1]*grid.cells[2],"Wrong result size");
            double mass=0;const double cellV=grid.spacing[0]*grid.spacing[1]*grid.spacing[2];
            for(vtkIdType i=0;i<rho->GetNumberOfTuples();++i){const double r=rho->GetComponent(i,0);check(std::isfinite(r)&&r>=0&&r<=1,"Invalid density");mass+=r*cellV;}
            const double error=std::abs(mass-reference)/reference;
            std::cout<<shape<<' '<<name<<" relative material error="<<error<<" time="<<seconds<<" s\n";
            const double dcSeconds=exportDC(amr,out/(shape+"_"+name+".obj"));
            timings<<shape<<" / "<<name<<" - "<<dcSeconds<<" s\n";
            timings.flush();check(bool(timings),"DC timings write failed");
            check(source->GetNumberOfCells()==count&&source->GetNumberOfPoints()==pointCount,"Source size changed");
            for(vtkIdType i=0;i<count;++i)check(source->GetCellData()->GetArray("vf")->GetComponent(i,0)==1.,"Source VF changed");
        }
    }
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}





