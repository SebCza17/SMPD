#include "mainwindow.h"
#include "ui_mainwindow.h"



#include <QImage>
#include <QDebug>
#include <boost/numeric/ublas/matrix.hpp>
#include "matrixutil.hpp"

using namespace boost::numeric::ublas;

#define myqDebug() qDebug() << string


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    FSupdateButtonState();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::updateDatabaseInfo()
{
    ui->FScomboBox->clear();
    for(unsigned int i=1; i<=database.getNoFeatures(); ++i)
        ui->FScomboBox->addItem(QString::number(i));

    ui->FStextBrowserDatabaseInfo->setText("noClass: " +  QString::number(database.getNoClass()));
    ui->FStextBrowserDatabaseInfo->append("noObjects: "  +  QString::number(database.getNoObjects()));
    ui->FStextBrowserDatabaseInfo->append("noFeatures: "  +  QString::number(database.getNoFeatures()));

}

void MainWindow::FSupdateButtonState(void)
{
    if(database.getNoObjects()==0)
    {
        FSsetButtonState(false);
    }
    else
        FSsetButtonState(true);

}


void MainWindow::FSsetButtonState(bool state)
{
   ui->FScomboBox->setEnabled(state);
   ui->FSpushButtonCompute->setEnabled(state);
   ui->FSpushButtonSaveFile->setEnabled(state);
   ui->FSradioButtonFisher->setEnabled(state);
   ui->FSradioButtonSFS->setEnabled(state);
}

void MainWindow::on_FSpushButtonOpenFile_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this,
        tr("Open TextFile"), "", tr("Texts Files (*.txt)"));

    if ( !database.load(fileName.toStdString()) )
        QMessageBox::warning(this, "Warning", "File corrupted !!!");
    else
        QMessageBox::information(this, fileName, "File loaded !!!");

    FSupdateButtonState();
    updateDatabaseInfo();
}

void MainWindow::on_FSpushButtonCompute_clicked()
{
    int dimension = ui->FScomboBox->currentText().toInt();


    if( ui->FSradioButtonFisher->isChecked())
    {
        if (dimension == 1 && database.getNoClass() == 2)
        {
            float FLD = 0, tmp;
            int max_ind = -1;

            for (uint i = 0; i < database.getNoFeatures(); ++i)
            {
                std::map<std::string, float> classAverages;

                std::map<std::string, float> classStds;

                for (auto const &ob : database.getObjects())
                {
                    classAverages[ob.getClassName()] += ob.getFeatures()[i];
                    classStds[ob.getClassName()] += ob.getFeatures()[i] * ob.getFeatures()[i];
                }

                std::for_each(database.getClassCounters().begin(), database.getClassCounters().end(), [&](const std::pair<std::string, int> &it)
                {
                    classAverages[it.first] /= it.second;

                    classStds[it.first] = std::sqrt(classStds[it.first] / it.second - classAverages[it.first] * classAverages[it.first]);}
                );
                tmp = std::abs(classAverages[database.getClassNames()[0]] - classAverages[database.getClassNames()[1]]) / (classStds[database.getClassNames()[0]] + classStds[database.getClassNames()[1]]);

                if (tmp > FLD)
                {
                    FLD = tmp;
                    max_ind = i;
                }
              }
            ui->FStextBrowserDatabaseInfo->append("max_ind: "  +  QString::number(max_ind) + " " + QString::number(FLD));
        }
        else if (dimension > 1 && database.getNoClass() == 2)
        {
            int d = dimension;
            int n = database.getNoFeatures();
            int* vectorF = new int[n];

            for(int i = 0; i < n; i++)
                vectorF[i] = i;

            std::map<int, std::map<std::string, double>> classAveragesTab;

            for (uint i = 0; i < database.getNoFeatures(); ++i){

                for (auto const &ob : database.getObjects()){
                   classAveragesTab[i][ob.getClassName()] += ob.getFeatures()[i];
                }

                std::for_each(database.getClassCounters().begin(), database.getClassCounters().end(), [&](const std::pair<std::string, int> &it){
                    classAveragesTab[i][it.first] /= it.second;
                });
            }

            printCombination(vectorF, n, d, classAveragesTab);

        }
    }
    else
    {
        if (database.getNoClass() == 2)
        {
            float FLD = 0, tmp;
            int top = -1;
            int* bestTab = new int[dimension];
            std::map<int, std::map<std::string, double>> classAveragesTab;
            int* bestData = new int[dimension];
            double max = 0;

            for(int d = 1; d <= dimension; d ++)
            {
                if(d == 1)
                {
                    for (uint i = 0; i < database.getNoFeatures(); ++i)
                    {
                        std::map<std::string, float> classStds;

                         for (auto const &ob : database.getObjects())
                         {
                             classAveragesTab[i][ob.getClassName()] += ob.getFeatures()[i];
                             classStds[ob.getClassName()] += ob.getFeatures()[i] * ob.getFeatures()[i];
                         }

                         std::for_each(database.getClassCounters().begin(), database.getClassCounters().end(), [&](const std::pair<std::string, int> &it)
                         {
                             classAveragesTab[i][it.first] /= it.second;

                             classStds[it.first] = std::sqrt(classStds[it.first] / it.second - classAveragesTab[i][it.first] * classAveragesTab[i][it.first]);
                         });

                         tmp = std::abs(classAveragesTab[i][database.getClassNames()[0]] - classAveragesTab[i][database.getClassNames()[1]]) / (classStds[database.getClassNames()[0]] + classStds[database.getClassNames()[1]]);

                         if (tmp > FLD)
                         {
                             FLD = tmp;
                             top = i;
                         }
                    }
                    bestTab[0] = top;
                    bestData[0] = top;
                }
                else
                {
                    int bestI = 0;
                    bestTab = bestData;
                    for(int i = 0; i < database.getNoFeatures(); i++)
                    {
                        bool isIn = false;

                        for(int in = 0; in < d; in++)
                            if(bestData[in] == i)
                                isIn = true;
                        if(!isIn){
                            bestTab[d - 1] = i;
                            double fisher = fisherFun(bestTab, d, classAveragesTab, max, bestData);

                            if(max < fisher){
                                max = fisher;
                                bestI = i;
                            }

                        }
                    }
                        bestData[d - 1] = bestI;
                }

            }




           qInfo() << max;

            for(int i = 0; i < dimension; i ++)
              qInfo() << bestData[i];

            qInfo() << "######";
       }
    }
}

void MainWindow::combinationUtil(int* arr, int* data, int start, int end, int index, int d, std::map<int, std::map<std::string, double>> classAveragesTab, double& max, int* bestData)  {
    if (index == d){

        double fisher = fisherFun(data, d, classAveragesTab, max, bestData);

        if(max < fisher){
            max = fisher;
            for(int i = 0; i < d; i++)
                bestData[i] = data[i];
        }

        return;
    }

    for (int i = start; i <= end && end - i + 1 >= d - index; i++){

        data[index] = arr[i];
        combinationUtil(arr, data, i+1, end, index+1, d, classAveragesTab, max, bestData);

    }

}

void MainWindow::printCombination(int* arr, int n, int d, std::map<int, std::map<std::string, double>> classAveragesTab){
    int* data = new int[d];
    double max = 0;
    int* bestData = new int[d];

    combinationUtil(arr, data, 0, n-1, 0, d, classAveragesTab, max, bestData);


    qInfo() << max;

    for(int i = 0; i < d; i ++)
        qInfo() << bestData[i];
}


double MainWindow::fisherFun(int* data, int d, std::map<int, std::map<std::string, double>> classAveragesTab, double& max, int* bestData)  {
    matrix<double> maS (d, 0);
    matrix<double> mbS (d, 0);

    matrix<double> maC (d, d);
    matrix<double> mbC (d, d);

    vector<double> vaU (d);
    vector<double> vbU (d);


    int j = 0, k = 0;

    for(int i = 0; i < d; i++){

        j = 0;
        k = 0;

        for (auto const &ob : database.getObjects()){

            if(ob.getClassName() == "Acer"){
                vaU(i) = classAveragesTab[data[i]][ob.getClassName()];

                if(maS.size2() <= j)
                    maS.resize(d, maS.size2() + 1);

                maS(i, j++) = ob.getFeatures()[data[i]] - vaU(i);


            }
            else{
                vbU(i) = classAveragesTab[data[i]][ob.getClassName()];

                if(mbS.size2() <= k)
                    mbS.resize(d, mbS.size2() + 1);
                mbS(i, k++) = ob.getFeatures()[data[i]] - vbU(i);
            }
        }


    }
    maC = prod(maS, trans(maS)) / j;
    mbC = prod(mbS, trans(mbS)) / k;

    //maC = prod(maS, trans(maS));
    //mbC = prod(mbS, trans(mbS));

    double detA = determinant(maC);
    double detB = determinant(mbC);

    vector<double> subU = (vaU - vbU);




    return sqrt(inner_prod(subU, subU)) / (detA + detB);;
}




void MainWindow::on_FSpushButtonSaveFile_clicked()
{
    QString fileName = QFileDialog::getSaveFileName(this,
    tr("Open TextFile"), "D:\\Users\\Krzysiu\\Documents\\Visual Studio 2015\\Projects\\SMPD\\SMPD\\Debug\\", tr("Texts Files (*.txt)"));

        QMessageBox::information(this, "My File", fileName);
        database.save(fileName.toStdString());
}

void MainWindow::on_PpushButtonSelectFolder_clicked()
{
}

void MainWindow::on_CpushButtonOpenFile_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this,
        tr("Open TextFile"), "", tr("Texts Files (*.txt)"));

    if ( !database.load(fileName.toStdString()) )
        QMessageBox::warning(this, "Warning", "File corrupted !!!");
    else
        QMessageBox::information(this, fileName, "File loaded !!!");

    FSupdateButtonStateC();
    updateDatabaseInfoC();
}

void MainWindow::on_CpushButtonSaveFile_clicked()
{

}

void MainWindow::on_CpushButtonTrain_clicked()
{
    int noTrainingPart = ui->CplainTextEditTrainingPart->toPlainText().toInt();

    if(noTrainingPart > 0 && noTrainingPart < database.getNoObjects() - 1 && ui && !ui->CheckBoxBoot->isChecked())
    {
        traningPartTab = new int[noTrainingPart];
        testPartTab = new int[database.getNoObjects() - noTrainingPart];

        for (int i = 0; i < noTrainingPart; i++)
        {
            traningPartTab[i] = rand() % (database.getNoObjects() - 1) + 1;
            for (int j = 0; j<i; j++)
            {
                if (traningPartTab[j] == traningPartTab[i])
                {
                    traningPartTab[i] = rand() % (database.getNoObjects() - 1) + 1;
                    j = -1;
                }
            }
        }

        int k = 0;
        for(int i = 0; i < database.getNoObjects(); i ++)
        {
            bool flag = true;
            for(int j = 0; j < noTrainingPart; j++)
            {
                if(traningPartTab[j] == i)
                    flag = false;
            }

            if(flag)
                testPartTab[k++] = i;
        }
    }
    else if (ui->CheckBoxBoot->isChecked())
    {
        int noN = ui->ComboBoxN->currentText().toInt();

        int* tmpTab = new int[database.getNoObjects()];

        nTrainingPartTab = new int*[noN];
        nTestPartTab = new int*[noN];
        sizeTab = new int[noN];

        for(int n = 0; n < noN; n ++)
        {
            int k = 0;
            std::map<int, int> tmpMap;
            for (int i = 0; i < database.getNoObjects(); i++)
            {
                tmpTab[i] = rand() % (database.getNoObjects() - 1) + 1;

                bool flag = true;

                for (int j = 0; j<i; j++)
                {
                    if(tmpTab[j] == tmpTab[i])
                        flag = false;
                    //qInfo() << tmpTab[j] << tmpTab[i];
                }
                if(flag){
                    tmpMap[k++] = tmpTab[i];

                }
            }

            nTrainingPartTab[n] = new int[k];
            for(int i = 0; i < k; i ++)
                nTrainingPartTab[n][i] = tmpMap[i];

            nTestPartTab[n] = new int[database.getNoObjects() - k];

            int p = 0;
            for(int i = 0; i < database.getNoObjects(); i ++)
            {
                bool flag = true;
                for(int j = 0; j < k; j++)
                {
                    if(nTrainingPartTab[n][j] == i)
                        flag = false;
                }

                if(flag)
                    nTestPartTab[n][p++] = i;
            }

            sizeTab[n] = k;
        }
    }

    ui->CpushButtonExecute->setEnabled(true);

    FSsetButtonStateC(false);
}
int MainWindow::methodNN(int noTrainingPart, int* testPartTab, int* trainingPartTab){
    int good = 0;

    for(int i = 0; i < database.getNoObjects() - noTrainingPart; i ++)
    {
        std::map<int, float> classObjectInner;
        int minKey = 0;

        for(int j = 0; j < noTrainingPart; j ++)
        {
            for(int k = 0; k < database.getNoFeatures(); k ++)
                classObjectInner[j] += (database.getObjects()[testPartTab[i]].getFeatures()[k] - database.getObjects()[trainingPartTab[j]].getFeatures()[k]) * (database.getObjects()[testPartTab[i]].getFeatures()[k] - database.getObjects()[trainingPartTab[j]].getFeatures()[k]);

            if(classObjectInner[minKey] > classObjectInner[j])
                minKey = j;
        }
            if(database.getObjects()[testPartTab[i]].getClassName() == database.getObjects()[traningPartTab[minKey]].getClassName())
                good ++;
    }

    return (good * 100) / (database.getNoObjects() - noTrainingPart);
}

int MainWindow::methodKNN(int noTrainingPart, int noK, int* testPartTab, int* trainingPartTab){
    int good = 0;

    for(int i = 0; i < database.getNoObjects() - noTrainingPart; i ++)
    {
        std::map<int, float> classObjectInner;

        int* minKeyTab = new int[noK];


        for(int j = 0; j < noTrainingPart; j ++)
        {
            for(int k = 0; k < database.getNoFeatures(); k ++)
                classObjectInner[j] += (database.getObjects()[testPartTab[i]].getFeatures()[k] - database.getObjects()[trainingPartTab[j]].getFeatures()[k]) * (database.getObjects()[testPartTab[i]].getFeatures()[k] - database.getObjects()[trainingPartTab[j]].getFeatures()[k]);

        }

        std::vector<std::pair<int, float>> pairs;
        for (auto itr = classObjectInner.begin(); itr != classObjectInner.end(); ++itr)
            pairs.push_back(*itr);

        sort(pairs.begin(), pairs.end(), [=](std::pair<int, float>& a, std::pair<int, float>& b)
        {
            return a.second < b.second;
        });


        for(int j = 0; j < noK; j ++)
            minKeyTab[j] = pairs[j].first;


        int isGood = 0;

        for(int l = 0; l < noK; l ++)
            if(database.getObjects()[testPartTab[i]].getClassName() == database.getObjects()[traningPartTab[minKeyTab[l]]].getClassName())
                isGood ++;
            else
                isGood --;

        if(isGood > 0)
            good ++;
    }

    return (good * 100) / (database.getNoObjects() - noTrainingPart);
}

int MainWindow::methodNM(int noTrainingPart, int* testPartTab, int* trainingPartTab){
    int good = 0;
    int countAcer = 0;

    std::map<std::string, std::map<int, float>> classObjectAvg;

    for(int j = 0; j < noTrainingPart; j ++)
    {
        for(int k = 0; k < database.getNoFeatures(); k ++)
            classObjectAvg[database.getObjects()[trainingPartTab[j]].getClassName()][k] += database.getObjects()[trainingPartTab[j]].getFeatures()[k];

        if(database.getObjects()[trainingPartTab[j]].getClassName() == "Acer")
            countAcer ++;
    }

    for(int k = 0; k < database.getNoFeatures(); k ++)
    {
        classObjectAvg["Acer"][k] /= countAcer;
        classObjectAvg["Quercus"][k] /= (noTrainingPart - countAcer);
    }


    for(int i = 0; i < database.getNoObjects() - noTrainingPart; i ++)
    {
        float minKey = 0;
        std::string minName;

        float divAcer;
        float divQuercus;

        for(int k = 0; k < database.getNoFeatures(); k ++)
        {
            divQuercus += (database.getObjects()[testPartTab[i]].getFeatures()[k] - classObjectAvg["Quercous"][k]) * (database.getObjects()[testPartTab[i]].getFeatures()[k] - classObjectAvg["Acer"][k]);
            divAcer += (database.getObjects()[testPartTab[i]].getFeatures()[k] - classObjectAvg["Acer"][k]) * (database.getObjects()[testPartTab[i]].getFeatures()[k] - classObjectAvg["Acer"][k]);

        }

        if(divAcer < divQuercus)
            minName = "Acer";
        else
            minName = "Quercus";

        if(database.getObjects()[testPartTab[i]].getClassName() == minName)
            good ++;
    }

    return (good * 100) / (database.getNoObjects() - noTrainingPart);
}

void MainWindow::on_CpushButtonExecute_clicked()
{
    int noN = ui->ComboBoxN->currentText().toInt();
    int noK = ui->CcomboBoxK->currentText().toInt();

    if(!ui->CheckBoxBoot->isChecked())
    {
        int noTrainingPart = ui->CplainTextEditTrainingPart->toPlainText().toInt();

        if(ui->CcomboBoxClassifiers->currentText() == "NN")
        {

            int percent = methodNN(noTrainingPart, testPartTab, traningPartTab);

            ui->CtextBrowser->append("NN Good: "  +  QString::number(percent) + "%");

        }
        else if(ui->CcomboBoxClassifiers->currentText() == "K-NN")
        {
            int percent = methodKNN(noTrainingPart, noK, testPartTab, traningPartTab);

            ui->CtextBrowser->append("K-NN Good: "  +  QString::number(percent) + "%");
        }
        else if(ui->CcomboBoxClassifiers->currentText() == "NM")
        {
            int percent = methodNM(noTrainingPart, testPartTab, traningPartTab);

            ui->CtextBrowser->append("NN Good: "  +  QString::number(percent) + "%");
        }
    }
    else if (ui->CheckBoxBoot->isChecked())
    {
        int percentAVG = 0;
        for(int n = 0; n < noN; n++)
        {

            int noTrainingPart  = sizeTab[n];

            int percent = 0;

            testPartTab = nTestPartTab[n];
            traningPartTab = nTrainingPartTab[n];

            if(ui->CcomboBoxClassifiers->currentText() == "NN")
            {

                percent = methodNN(noTrainingPart, testPartTab, traningPartTab);

                ui->CtextBrowser->append("NN Good: "  +  QString::number(percent) + "%");

            }
            else if(ui->CcomboBoxClassifiers->currentText() == "K-NN")
            {
                percent = methodKNN(noTrainingPart, noK, testPartTab, traningPartTab);

                ui->CtextBrowser->append("K-NN Good: "  +  QString::number(percent) + "%");
            }
            else if(ui->CcomboBoxClassifiers->currentText() == "NM")
            {

                percent =  methodNM(noTrainingPart, testPartTab, traningPartTab);

                ui->CtextBrowser->append("NN Good: "  +  QString::number(percent) + "%");
            }

            percentAVG += percent;
        }

        ui->CtextBrowser->append("NN Good: N-times"  +  QString::number(percentAVG/noN) + "%");
    }



}

void MainWindow::updateDatabaseInfoC()
{
    ui->CcomboBoxClassifiers->clear();
    ui->CcomboBoxClassifiers->addItem("NN");
    ui->CcomboBoxClassifiers->addItem("K-NN");
    ui->CcomboBoxClassifiers->addItem("NM");
    ui->CcomboBoxClassifiers->addItem("K-NN");



    ui->ComboBoxN->clear();
    for(unsigned int i=1; i<= database.getNoFeatures(); ++i)
        ui->ComboBoxN->addItem(QString::number(i));


    ui->CcomboBoxK->clear();
    for(unsigned int i=1; i<=database.getNoFeatures(); i += 2)
        ui->CcomboBoxK->addItem(QString::number(i));

    ui->CtextBrowser->setText("noClass: " +  QString::number(database.getNoClass()));
    ui->CtextBrowser->append("noObjects: "  +  QString::number(database.getNoObjects()));
    ui->CtextBrowser->append("noFeatures: "  +  QString::number(database.getNoFeatures()));

}

void MainWindow::FSupdateButtonStateC(void)
{
    if(database.getNoObjects()==0)
    {
        FSsetButtonStateC(false);
    }
    else
        FSsetButtonStateC(true);

}


void MainWindow::FSsetButtonStateC(bool state)
{
   ui->CpushButtonTrain->setEnabled(state);
   ui->CplainTextEditTrainingPart->setEnabled(state);

}
