#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

QT_FORWARD_DECLARE_CLASS(Generator)

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	explicit MainWindow(QWidget *parent = 0);
	~MainWindow();
private slots:
	void updateFieldOfViewDisplay();
private:
	Ui::MainWindow *ui;
	Generator * m_generator;
};

#endif // MAINWINDOW_H
