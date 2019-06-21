#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <fstream>
#include <QCheckBox>
#include <QFileDialog>
#include <QtCore>
#include <QtGui>
#include <QMessageBox>


#include"database.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    void printCombination(int* arr, int n, int d, std::map<int, std::map<std::string, double>>);
    void combinationUtil(int* arr, int* data, int start, int end, int index, int d, std::map<int, std::map<std::string, double>> classAveragesTab, double& max, int* bestData);
    double fisherFun(int* data, int d, std::map<int, std::map<std::string, double>> classAveragesTab, double& max, int* bestData);
    int methodNN(int noTrainingPart, int* testPartTab, int* trainingPartTab);
    int methodKNN(int noTrainingPart, int noK, int* testPartTab, int* trainingPartTab);
    int methodNM(int noTrainingPart, int* testPartTab, int* trainingPartTab);

private:
    bool loadFile(const std::string &fileName);
    void updateDatabaseInfo();
    void saveFile(const std::string &fileName);

    void FSupdateButtonState(void);
    void FSsetButtonState(bool state);

    void updateDatabaseInfoC();
    void FSupdateButtonStateC(void);
    void FSsetButtonStateC(bool state);




private slots:
    void on_FSpushButtonOpenFile_clicked();

    void on_FSpushButtonCompute_clicked();

    void on_FSpushButtonSaveFile_clicked();

    void on_PpushButtonSelectFolder_clicked();


    void on_CpushButtonOpenFile_clicked();

    void on_CpushButtonSaveFile_clicked();

    void on_CpushButtonTrain_clicked();

    void on_CpushButtonExecute_clicked();

private:
    Ui::MainWindow *ui;

private:
     Database database;
     int* traningPartTab;
     int* testPartTab;
     int** nTrainingPartTab;
     int** nTestPartTab;
     int* sizeTab;
};

#endif // MAINWINDOW_H
