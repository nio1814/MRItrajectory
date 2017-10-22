#ifndef PHANTOMRECONSTRUCTION_H
#define PHANTOMRECONSTRUCTION_H

#include <QWidget>
#include <QPointer>

#include <memory>

QT_FORWARD_DECLARE_CLASS(Phantom)
QT_FORWARD_DECLARE_STRUCT(Trajectory)
QT_FORWARD_DECLARE_CLASS(QLabel)

class PhantomReconstruction : public QWidget
{
	Q_OBJECT
public:
	explicit PhantomReconstruction(QWidget *parent = 0);

signals:

public slots:
	void reconstruct(Trajectory *trajectory);
protected:
	void paintEvent(QPaintEvent *event);
private:
	std::shared_ptr<Phantom> m_phantom2D;
	QPointer<QLabel> m_imageLabel;
};

#endif // PHANTOMRECONSTRUCTION_H
